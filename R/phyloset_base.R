
quantile_rank <- function(x) {
    ranks <- base::rank(x, ties.method = "average")
    (ranks - 0.5) / length(x)
}

TI_map <- list(TXI = "Transcriptomic Index",
               TAI = "Transcriptomic Age Index",
               TDI = "Transcriptomic Divergence Index",
               TPI = "Transcriptomic Polymorphism Index",
               TEI = "Transcriptomic Evolutionary Index"
               )


#' @import S7
PhyloExpressionSet <- new_class("PhyloExpressionSet",
    properties = list(
        ## CONSTRUCTOR PARAMETERS
        data = new_required_property( # the full phyex set 
            class = class_data.frame,
            name = "data"
            ),
        index_type = new_options_property( # the type of transcriptomic index
            class = class_character,
            options = names(TI_map),
            default = "TXI"
            ),
        conditions_label = new_property(
            class = class_character,
            default = "Ontogeny"
            ),
        bootstrap_sample_size = new_property(
            class = class_numeric,
            default = 5000L
            ),
        null_conservation_sample_size = new_property(
            class = class_numeric,
            default = 5000L
            ),
        ## FIELDS & PROPERTIES
        index_full_name = new_property(
            class = class_character,
            getter <- \(self) TI_map[[self@index_type]]
            ),
        strata_vector = new_property(
            class = class_double,
            getter = function(self) { 
                x <- self@data[[1]]
                return(x)
            },
            setter = function(self, value) {
                if (!length(value))
                    return(self)
                self@data[[1]] <- value
                self
            }
            ),
        gene_ids = new_property(
            class = class_factor,
            getter = \(self) self@data[[2]],
            setter = function(self, value) {
                if (!length(value))
                    return(self)
                self@data[[2]] <- value
                self
            }
            ),
        count_matrix = new_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            getter = function(self) {
                m <- self@data[3:ncol(self@data)] |>
                    as.matrix()
                rownames(m) <- self@gene_ids
                
                return(m)
            },
            setter = function(self, value) {
                if (!length(value))
                    return(self)
                self@data[3:ncol(self@data)] <- value
                self
            },
            validator = \(value) if (any(is.na(value))) "cannot contain NA values. Check data[3:ncol(data)]."
            ),
        conditions = new_property(
            class = class_factor,
            getter = \(self) factor(colnames(self@count_matrix), 
                                    levels=unique(colnames(self@count_matrix),
                                    ordered=TRUE))
            ),
        num_genes = new_property(
            class = class_integer,
            getter = \(self) nrow(self@count_matrix)
            ),
        num_conditions = new_property(
            class = class_integer,
            getter = \(self) ncol(self@count_matrix)
            ),
        pTXI = new_cached_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            getter = function(self) {
                m <- cpp_pMatrix(self@count_matrix, self@strata_vector)
                rownames(m) <- self@gene_ids
                colnames(m) <- self@conditions
                return(m)
            }
            ),
        TXI = new_property(
            class = class_double,
            getter = \(self) colSums(self@pTXI)
            ),
        bootstrapped_txis = new_cached_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            getter = .generate_bootstrapped_txis
            ),
        null_conservation_txis = new_cached_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            getter = .generate_conservation_txis
        )
        )
    
    )

# S7::method(print, PhyloExpressionSet) <- function(x) {
#     return(print(x@TXI))
# }



TXI_conf_int <- function(phyex_set, 
                         probs=c(.025, .975)) {
    CIs <- apply(phyex_set@bootstrapped_txis, 2, quantile, probs=probs)
    return (CIs)
}

sTXI <- function(phyex_set,
                 option="identity") {
    mat <- rowsum(phyex_set@pTXI, phyex_set@strata_vector)
    if (option == "add") {
        mat <- apply(mat,2,cumsum)
    }
    if (! option %in% c("identity", "add"))
        stop("'option' must be one of 'identity' or 'add",
             call. = FALSE
        )
    return(mat)
}

.generate_bootstrapped_txis <- function(phyex_set,
                                        sample_size=NULL) {
    if (is.null(sample_size))
        sample_size <- phyex_set@null_conservation_sample_size
    
    N <- phyex_set@num_genes
    
    sampled_indices <- matrix(sample(N, size=N*sample_size, replace=TRUE), 
                              nrow=N, ncol=sample_size)
    selection_matrix <- Matrix::sparseMatrix(i = as.vector(sampled_indices),
                                             j = rep(1:sample_size, each = N),
                                             x = 1,
                                             dims= c(N, sample_size))
    
    bootstrap_matrix <- t(phyex_set@pTXI) %*% selection_matrix |> 
        t() # columns: conditions. rows: bootstraps
    
    colnames(bootstrap_matrix) <- phyex_set@conditions
    
    return(bootstrap_matrix)
}

.generate_conservation_txis <- function(phyex_set, sample_size=NULL) {
    if (is.null(sample_size))
        sample_size <- phyex_set@null_conservation_sample_size
    permuted_txi_matrix <- cpp_bootMatrix(phyex_set@count_matrix,
                                          phyex_set@strata_vector,
                                          sample_size)
    
    colnames(permuted_txi_matrix) <- phyex_set@conditions
    return(permuted_txi_matrix)
}


#' #' @import S7
#' PerturbedPhyloExpressionSet <- new_class("PerturbedPhyloExpressionSet",
#'     parent = PhyloExpressionSet,
#'     properties = list(
#'         removed_genes = new_property(
#'             class = class_character,
#'             default = character(0)
#'             ),
#'         genes_bitvector = new_property(
#'             class = class_logical,
#'             getter = \(self) ! self@gene_ids %in% self@removed_genes
#'             ),
#'         TXI = new_property(
#'             class = class_double,
#'             getter = \(self) matrixStats::colSums2(self@pTXI, rows = which(self@genes_bitvector)) # using colSums2 avoids unnecessary copying
#'             )
#'         )
#'     )

