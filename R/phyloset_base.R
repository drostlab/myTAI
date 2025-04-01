
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
        ## PARAMETERS
        #TODO add species and process/ title of data
        strata_vector = new_required_property(
            class = class_numeric,
            validator = \(value) if (any(is.na(value))) "cannot contain NA values. Check data[1].",
            name = "strata_vector"
        ),
        gene_ids = new_required_property(
            class = class_character,
            validator = \(value) if (any(is.na(value))) "cannot contain NA values. Check data[2].",
            name = "gene_ids"
        ),
        count_matrix = new_required_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            validator = \(value) if (any(is.na(value))) "cannot contain NA values. Check data[3:ncol(data)].",
            name = "count_matrix"
        ),
        name = new_property(
            class = class_character,
            default = "Phylo Expression Set"
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
        is_time_series = new_property(
            class = class_logical,
            default = TRUE
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
        data = new_property(
            class = class_data.frame,
            getter <- \(self) tibble::tibble(Stratum=self@strata_vector, 
                                             GeneID=self@gene_ids, 
                                             tibble::as_tibble(self@count_matrix))
        ),
        index_full_name = new_property(
            class = class_character,
            getter <- \(self) TI_map[[self@index_type]]
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
        bootstrapped_txis = new_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            getter = \(self) memo_generate_bootstrapped_txis(self@pTXI,
                                                             self@count_matrix,
                                                             self@bootstrap_sample_size)
            ),
        null_conservation_txis = new_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            getter = \(self) memo_generate_conservation_txis(self@strata_vector,
                                                             self@count_matrix,
                                                             self@null_conservation_sample_size)
        )
        ),
    constructor = function (data = NULL, 
                            strata_vector = as.numeric(data[[1]]),
                            gene_ids = as.character(data[[2]]),
                            count_matrix = as.matrix(data[3:ncol(data)]),
                            name = deparse(substitute(data)),
                            index_type = "TXI", 
                            conditions_label = "Ontogeny", 
                            is_time_series = TRUE,
                            bootstrap_sample_size = 5000L, 
                            null_conservation_sample_size = 5000L) {

        new_object(S7_object(), 
                   strata_vector = strata_vector,
                   gene_ids = gene_ids,
                   count_matrix = count_matrix,
                   name = name,
                   index_type = index_type, 
                   conditions_label = conditions_label, 
                   is_time_series = is_time_series,
                   bootstrap_sample_size = bootstrap_sample_size, 
                   null_conservation_sample_size = null_conservation_sample_size)
        }
    
    )


S7::method(print, PhyloExpressionSet) <- function(x, ...) {
    cat("\n", "Phylo Expression Set", "\n", sep="")
    cat(x@conditions_label, ":  ", paste(as.character(x@conditions), collapse = ", "), "\n", sep="")
    cat("Number of genes: ", x@num_genes, "\n", sep="")
}



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

transform_counts <- S7::new_generic("transform_counts", "phyex_set")
S7::method(transform_counts, PhyloExpressionSet) <- function(phyex_set, 
                                                             FUN,
                                                             FUN_name=deparse(substitute(FUN)),
                                                             new_name=paste(phyex_set@name, "transformed by", FUN_name)) {
    f <- match.fun(FUN)
    phyex_set@count_matrix <- f(phyex_set@count_matrix)
    phyex_set@name <- new_name
    return(phyex_set)
}

select_genes <- S7::new_generic("select_genes", "phyex_set")
S7::method(select_genes, PhyloExpressionSet) <- function(phyex_set, 
                                                         genes) {
    indices <- (phyex_set@gene_ids %in% genes)
    
    phyex_set@strata_vector <- phyex_set@strata_vector[indices]
    phyex_set@gene_ids <- phyex_set@gene_ids[indices]
    phyex_set@count_matrix <- phyex_set@count_matrix[indices, ]

    return(phyex_set)
}

remove_genes <- function(phyex_set, genes) {
    selected_genes <- setdiff(phyex_set@gene_ids, genes)
    select_genes(phyex_set, selected_genes)
}







