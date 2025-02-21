
quantile_rank <- function(x) {
    ranks <- base::rank(x, ties.method = "average")
    (ranks - 0.5) / length(x)
}

TI_map = list(TAI = "Transcriptomic Age Index",
                   TDI = "Transcriptomic Divergence Index",
                   TPI = "Transcriptomic Polymorphism Index",
                   TEI = "Transcriptomic Evolutionary Index",
                   TXI = "Transcriptomic Index")


#' @import S7
PhyloExpressionSet <- new_class("PhyloExpressionSet",
    properties = list(
        ## CONSTRUCTOR PARAMETERS
        data = new_property( # the full phyex set 
            class = class_data.frame,
            default = quote(stop("@data is required"))
            ),
        index_type = new_property( # the type of transcriptomic index
            class = class_character,
            default = "TXI",
            validator = \(value) if (! value %in% names(TI_map)) paste("must be one of", paste(names(TI_map), collapse=", "))
            ),
        stages_label = new_property(
            class = class_character,
            default = "Ontogeny"
            ),
        strata_transform = new_property(
            class = class_function,
            default = quantile_rank
            ),
        count_transform = new_property(
            class = class_function,
            default = sqrt
            ),
        bootstrap_sample_size = new_property(
            class = class_integer,
            default = 5000L
            ),
        null_conservation_sample_size = new_property(
            class = class_integer,
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
                x <- self@data[[1]] |>
                    self@strata_transform() # apply strata transformation
                return(x)
            }
            ),
        gene_ids = new_property(
            class = class_factor,
            getter = \(self) self@data[[2]]
            ),
        count_matrix = new_property(
            class = class_double,
            getter = function(self) {
                m <- self@data[3:ncol(self@data)] |>
                    as.matrix() |>
                    self@count_transform() # apply transformations
                rownames(m) <- self@gene_ids
                
                return(m)
            },
            validator = \(value) if (any(is.na(value))) "cannot contain NA values. Check data[3:ncol(data)]."
            ),
        stages = new_property(
            class = class_factor,
            getter = \(self) factor(colnames(self@count_matrix), 
                                    levels=unique(colnames(self@count_matrix),
                                    ordered=TRUE))
            ),
        num_genes = new_property(
            class = class_integer,
            getter = \(self) nrow(self@count_matrix)
            ),
        num_stages = new_property(
            class = class_integer,
            getter = \(self) ncol(self@count_matrix)
            ),
        pTXI = new_property(
            class = class_double,
            getter = function(self) {
                m <- cpp_pMatrix(self@count_matrix, self@strata_vector)
                rownames(m) <- self@gene_ids
                colnames(m) <- self@stages
                return(m)
            }
            ),
        TXI = new_property(
            class = class_double,
            getter = \(self) colSums(self@pTXI)
            ),
        bootstrapped_TXIs = new_property(
            class = class_double,
            getter = \(self) generate_bootstrapped_txi(self@pTXI, self@bootstrap_sample_size)
            ),
        null_conservation_TXIs = new_property(
            class = class_double,
            getter = function(self) {
                permuted_txi_matrix <- cpp_bootMatrix(self@count_matrix,
                                                      self@strata_vector,
                                                      self@null_conservation_sample_size)
                
                colnames(permuted_txi_matrix) <- self@stages
                return(permuted_txi_matrix)
            }
        )
        )
    
    )

# PLOTTING TAI

TXI_conf_int <- function(phyex_set, 
                         probs=c(.025, .975)) {
    CIs <- apply(phyex_set@bootstrapped_TXIs, 2, quantile, probs=probs)
    return (CIs)
}

generate_bootstrapped_txi <- function(pTXI,
                                      num_bootstraps=2000) {
    num_genes <- nrow(pTXI)
    freq_matrix <- replicate(num_bootstraps, sample(num_genes, replace = TRUE)) |> # sample genes with replacement
        apply(2, \(col) tabulate(col, nbins=num_genes)) # convert matrix to frequency matrix
        
    bootstrap_matrix <- t(pTXI) %*% freq_matrix |> 
        t() # columns: stages. rows: bootstraps
    colnames(bootstrap_matrix) <- colnames(pTXI)
    
    return(bootstrap_matrix)
}



#' @import ggplot2 tibble
plot_bootstrap_sample <- function(phyex_set) {
    df_sample <- as_tibble(phyex_set@bootstrapped_TXIs) |>
        rowid_to_column("Id") |>
        tidyr::pivot_longer(cols=phyex_set@stages, names_to="Stage", values_to = "Value")
    df_true <- tibble::tibble(Stage = phyex_set@stages,
                              TXI = phyex_set@TXI
                              )
    
    p <- ggplot(data=df_sample) +
        geom_line(
            aes(
                x=factor(Stage, levels=unique(Stage)),
                y=Value,
                group=Id
            ),
            alpha=0.2,
            colour="gray70"
        ) +
        geom_line(
            data=df_true,
            aes(
                x=Stage,
                y=TXI,
                group=0),
            lwd=2,
            alpha=1.0,
            colour="black"
        ) +
        labs(
            x=phyex_set@stages_label,
            y=phyex_set@index_full_name
        ) +
        theme_minimal()
    
    return(p)
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


#' @import S7
PerturbedPhyloExpressionSet <- new_class("PerturbedPhyloExpressionSet",
    parent = PhyloExpressionSet,
    properties = list(
        removed_genes = new_property(
            class = class_character,
            default = character(0)
            ),
        genes_bitvector = new_property(
            class = class_logical,
            getter = \(self) ! self@gene_ids %in% self@removed_genes
            ),
        TXI = new_property(
            class = class_double,
            getter = \(self) matrixStats::colSums2(self@pTXI, rows = which(self@genes_bitvector)) # using colSums2 avoids unnecessary copying
            )
        )
    )

