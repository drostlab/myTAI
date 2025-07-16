#' @title Generate Phylostratum Colors
#' @description Generate a color palette for phylostrata visualization using a log-scaled transformation.
#' @param n number of colors to generate
#' @return A character vector of color codes
#' @examples
#' # Generate colors for 5 phylostrata
#' colors <- PS_colours(5)
#' @importFrom grDevices colorRampPalette
#' @export
PS_colours <- function(n) {
    vals <- 1:n |>
        log() |>
        scales::rescale()
   
    pal <- grDevices::colorRampPalette(c("black", "#AD6F3B", "lightgreen"))
   
    pal(100)[floor(vals * 99) +1]
}   



#' @title Transcriptomic Index Name Mapping
#' @description Named list mapping transcriptomic index abbreviations to full names.
#' @format A named list with 5 elements:
#' \describe{
#'   \item{TXI}{Transcriptomic Index}
#'   \item{TAI}{Transcriptomic Age Index}
#'   \item{TDI}{Transcriptomic Divergence Index}
#'   \item{TPI}{Transcriptomic Polymorphism Index}
#'   \item{TEI}{Transcriptomic Evolutionary Index}
#' }
TI_map <- list(TXI = "Transcriptomic Index",
               TAI = "Transcriptomic Age Index",
               TDI = "Transcriptomic Divergence Index",
               TPI = "Transcriptomic Polymorphism Index",
               TEI = "Transcriptomic Evolutionary Index")


#' @title PhyloExpressionSet S7 Class
#' @description S7 class for storing and manipulating phylotranscriptomic expression data.
#' This class integrates phylostratum information with gene expression data across 
#' developmental stages or conditions.
#' 
#' @slot strata Factor vector of phylostratum assignments for each gene
#' @slot gene_ids Character vector of gene identifiers
#' @slot counts Matrix of expression counts with genes as rows and samples as columns
#' @slot groups Factor vector indicating which group each sample belongs to
#' @slot name Character string naming the dataset (default: "Phylo Expression Set")
#' @slot species Character string indicating the species (default: NULL)
#' @slot index_type Character string indicating the type of transcriptomic index (default: "TXI")
#' @slot conditions_label Character string labeling the conditions (default: "Ontogeny")
#' @slot is_time_series Logical indicating if data represents a time series (default: TRUE)
#' @slot null_conservation_sample_size Numeric sample size for null conservation tests (default: 5000)
#' 
#' @details
#' The PhyloExpressionSet class provides computed properties including:
#' - `counts_collapsed`: Expression data collapsed across replicates
#' - `data`: Complete dataset as a tibble
#' - `pTXI`: Phylostratum-specific transcriptomic index values
#' - `TXI`: Overall transcriptomic index per condition
#' - `conditions`: Factor of condition names
#' - `num_genes`, `num_conditions`, `num_strata`: Counts of dataset dimensions
#' 
#' @examples
#' # Create a PhyloExpressionSet from data
#' # phyex_set <- as_PhyloExpressionSet(my_data)
#' 
#' @import S7
#' @export
PhyloExpressionSet <- new_class("PhyloExpressionSet",
    properties = list(
        ## PARAMETERS
        # REQUIRED
        strata = new_required_property(
            class = class_factor,
            validator = function(value) {
                if (any(is.na(value))) "cannot contain NA values. Check data[1]."
                if (length(value) == 0) "cannot be empty. Check data[1]"
            },
            name = "strata"
        ),
        gene_ids = new_required_property(
            class = class_character,
            validator = function(value) {
                if (any(is.na(value))) "cannot contain NA values. Check data[2]."
                if (length(value) == 0) "cannot be empty. Check data[2]"
            },
            name = "gene_ids"
        ),
        counts = new_required_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            validator = function(value) {
                if (any(is.na(value))) "cannot contain NA values. Check data[3:ncol(data)]."
                if (length(value) == 0) "cannot be empty. Check data[3:ncol(data)]"
            },
            name = "counts"
        ),
        groups = new_required_property(
            class = class_factor,
            name = "groups"
        ),
        # OPTIONAL
        name = new_property(
            class = class_character,
            default = "Phylo Expression Set"
        ),
        species = new_property(
            class = class_character,
            default = NULL
        ),
        index_type = new_options_property( # the type of transcriptomic indexclass = class_character,
            options = names(TI_map),
            default = "TXI"),
        conditions_label = new_property(
            class = class_character,
            default = "Ontogeny"
            ),
        is_time_series = new_property(
            class = class_logical,
            default = TRUE
            ),
        null_conservation_sample_size = new_property(
            class = class_numeric,
            default = 5000L
            ),
        ## FIELDS & PROPERTIES
        counts_collapsed = new_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            getter = \(self) .collapse_replicates(self@counts, self@groups),
            validator = \(value) if (any(is.na(value))) "cannot contain NA values. Check data[3:ncol(data)].",
        ),
        data = new_property(
            class = class_data.frame,
            getter <- \(self) tibble::tibble(Stratum=self@strata, 
                                             GeneID=self@gene_ids, 
                                             tibble::as_tibble(self@counts))
        ),
        data_collapsed = new_property(
            class = class_data.frame,
            getter <- \(self) tibble::tibble(Stratum=self@strata, 
                                             GeneID=self@gene_ids, 
                                             tibble::as_tibble(self@counts_collapsed))
        ),
        index_full_name = new_property(
            class = class_character,
            getter <- \(self) TI_map[[self@index_type]]
        ),
        conditions = new_property(
            class = class_factor,
            getter = \(self) factor(colnames(self@counts_collapsed), 
                                    levels=unique(colnames(self@counts_collapsed),
                                    ordered=TRUE))
        ),
        num_genes = new_property(
            class = class_integer,
            getter = \(self) nrow(self@counts)
        ),
        num_conditions = new_property(
            class = class_integer,
            getter = \(self) length(self@conditions)
        ),
        num_strata = new_property(
            class = class_integer,
            getter = \(self) length(unique(self@strata))
        ),
        pTXI = new_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            getter = function(self) {
                m <- pTXI(self@counts_collapsed, self@strata)
                rownames(m) <- self@gene_ids
                colnames(m) <- self@conditions
                return(m)
            }
        ),
        TXI = new_property(
            class = class_double,
            getter = \(self) colSums(self@pTXI)
        ),
        precomputed_null_conservation_txis = new_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            default = NULL
        ),
        null_conservation_txis = new_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            getter = function(self) {
                if (is.null(self@precomputed_null_conservation_txis))
                    memo_generate_conservation_txis(self@strata,
                                                   self@counts_collapsed,
                                                   self@null_conservation_sample_size)
                else
                    self@precomputed_null_conservation_txis
            }
        ),
        sample_names = new_property(
            class = class_character,
            getter = function(self) colnames(self@counts)
        ),
        group_map = new_property(
            class = class_list,
            getter = function(self) split(self@sample_names, self@groups)
        ),
        pTXI_reps = new_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            getter = function(self) {
                m <- pTXI(self@counts, self@strata)
                rownames(m) <- self@gene_ids
                colnames(m) <- self@sample_names
                return(m)
            }
        ),
        TXI_reps = new_property(
            class = class_double,
            getter = \(self) colSums(self@pTXI_reps)
        )
        )
    )


#' @title Convert Data to PhyloExpressionSet
#' @description Convert a data frame with phylostratum, gene ID, and expression data 
#' into a PhyloExpressionSet object.
#' 
#' @param data A data frame where column 1 contains phylostratum information, 
#' column 2 contains gene IDs, and columns 3+ contain expression data
#' @param groups A factor or character vector indicating which group each sample belongs to.
#' Default uses column names from expression data
#' @param name A character string naming the dataset. Default uses the variable name
#' @param strata_labels Optional character vector of labels for phylostrata. 
#' If NULL, uses sorted unique values from column 1
#' @param ... Additional arguments passed to PhyloExpressionSet constructor
#' 
#' @return A PhyloExpressionSet object
#' 
#' @examples
#' # Convert data frame to PhyloExpressionSet
#' # phyex_set <- as_PhyloExpressionSet(my_data, 
#' #                                   groups = c("stage1", "stage1", "stage2", "stage2"),
#' #                                   name = "Development Dataset")
#' 
#' @export
as_PhyloExpressionSet <- function(data, 
                                  groups = colnames(data[,3:ncol(data)]),
                                  name = deparse(substitute(data)),
                                  strata_labels = NULL,
                                  ...) {
    gene_ids <- as.character(data[[2]])
    if (is.null(strata_labels))
        strata_labels <- sort(unique(as.numeric(data[[1]])))
    strata <- factor(as.numeric(data[[1]]), levels=sort(unique(as.numeric(data[[1]]))), labels=strata_labels)
    names(strata) <- gene_ids

    groups <- factor(groups, levels=unique(groups))
    
    counts <- as.matrix(data[3:ncol(data)])
    rownames(counts) <- gene_ids
    
    return(PhyloExpressionSet(
        strata = strata,
        gene_ids = gene_ids,
        counts = counts,
        groups = groups,
        name = name,
        ...
    ))
}

#' @title Match Gene Expression Data with Phylostratum Map
#' @description Join gene expression data with a phylostratum mapping to create 
#' a PhyloExpressionSet object.
#' 
#' @param data A data frame where column 1 contains gene IDs and columns 2+ contain expression data
#' @param phylomap A data frame with two columns: phylostratum assignments and gene IDs
#' @param groups A factor or character vector indicating which group each sample belongs to.
#' Default uses column names from expression data
#' @param name A character string naming the dataset. Default uses the variable name
#' @param ... Additional arguments passed to as_PhyloExpressionSet
#' 
#' @return A PhyloExpressionSet object
#' 
#' @examples
#' # Match expression data with phylostratum map
#' # phyex_set <- match_map(expression_data, phylo_map, 
#' #                        groups = c("stage1", "stage2", "stage3"),
#' #                        name = "Matched Dataset")
#' 
#' @importFrom dplyr inner_join relocate
#' @export
match_map <- function(data, 
                      phylomap,
                      groups = colnames(data[,2:ncol(data)]),
                      name = NULL,
                      ...) {
    if (is.null(name)) name <- deparse(substitute(data))
    colnames(phylomap) <- c("Stratum", "GeneID")
    
    data <- data |>
        inner_join(phylomap, by="GeneID") |>
        relocate(Stratum, .before = 1)
    
    return(as_PhyloExpressionSet(data, groups = groups, name = name, ...))
}




S7::method(print, PhyloExpressionSet) <- function(x, ...) {
    cat("\n", "Phylo Expression Set", "\n", sep="")
    cat("\n", x@name, "\n", sep="")
    cat(x@conditions_label, ":  ", paste(as.character(x@conditions), collapse = ", "), "\n", sep="")
    cat("Number of genes: ", x@num_genes, "\n", sep="")
}


#' @title Collapse Expression Data Across Replicates
#' @description Internal function to collapse expression data across replicates by taking row means.
#' 
#' @param counts Matrix of expression counts with genes as rows and samples as columns
#' @param groups Factor or character vector indicating which group each sample belongs to
#' 
#' @return Matrix with collapsed expression data, one column per group
#' 
#' @keywords internal
.collapse_replicates <- function(counts, groups) {
    m <- do.call(cbind, lapply(unique(groups), \(g) rowMeans(counts[, groups == g, drop = FALSE])))
    colnames(m) <- unique(groups)
    return(m)
}

#' @title Collapse PhyloExpressionSet Replicates
#' @description Convert a PhyloExpressionSet with replicates to one with collapsed expression data.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' 
#' @return A new PhyloExpressionSet object with collapsed expression data
#' 
#' @examples
#' # Collapse replicates in a PhyloExpressionSet
#' # collapsed_set <- collapse(phyex_set)
#' 
#' @export
collapse <- function(phyex_set) {
    data <- tibble::tibble(Stratum=phyex_set@strata, 
                           GeneID=phyex_set@gene_ids, 
                           tibble::as_tibble(phyex_set@counts_collapsed))
    as_PhyloExpressionSet(data)
}

#' @title Normalize Stage Expression Data
#' @description Normalize expression data to a specified total expression level per sample.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param total Numeric value to normalize each sample to (default: 1e6)
#' 
#' @return A PhyloExpressionSet object with normalized expression data
#' 
#' @examples
#' # Normalize to 1 million total expression per sample
#' # normalized_set <- normalise_stage_expression(phyex_set, total = 1e6)
#' 
normalise_stage_expression <- function(phyex_set, total=1e6) {
    phyex_set@counts <- sweep(phyex_set@counts, 2, colSums(phyex_set@counts), FUN="/") * total
    phyex_set
}

#' @title Transform Expression Counts in PhyloExpressionSet
#' @description Apply a transformation function to the expression counts in a PhyloExpressionSet.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param FUN Function to apply to the counts matrix
#' @param FUN_name Character string naming the transformation function (default: derived from FUN)
#' @param new_name Character string for the new dataset name (default: auto-generated)
#' @param ... Additional arguments passed to the transformation function
#' 
#' @return A PhyloExpressionSet object with transformed expression data
#' 
#' @examples
#' # Apply log transformation
#' # log_set <- transform_counts(phyex_set, log1p, "log1p")
#' 
#' @export
transform_counts <- S7::new_generic("transform_counts", "phyex_set")

#' @export
S7::method(transform_counts, PhyloExpressionSet) <- function(phyex_set, 
                                                             FUN,
                                                             FUN_name=deparse(substitute(FUN)),
                                                             new_name=paste(phyex_set@name, "transformed by", FUN_name),
                                                             ...) {
    f <- match.fun(FUN)
    phyex_set@counts <- f(phyex_set@counts, ...)
    phyex_set@name <- new_name
    return(phyex_set)
}

#' @title Short Alias for Transform Counts
#' @description Convenience alias for transform_counts function.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param ... Arguments passed to transform_counts
#' 
#' @return A PhyloExpressionSet object with transformed expression data
#' 
#' @examples
#' # Short alias for transformation
#' # log_set <- tf(phyex_set, log1p)
#' 
#' @export
tf <- transform_counts

#' @title Select Genes from PhyloExpressionSet
#' @description Extract a subset of genes from a PhyloExpressionSet object.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param genes Character vector of gene IDs to select
#' @param ... Additional arguments (currently unused)
#' 
#' @return A PhyloExpressionSet object containing only the selected genes
#' 
#' @examples
#' # Select specific genes
#' # selected_set <- select_genes(phyex_set, c("gene1", "gene2", "gene3"))
#' 
#' @export
select_genes <- S7::new_generic("select_genes", "phyex_set")

#' @export
S7::method(select_genes, PhyloExpressionSet) <- function(phyex_set, 
                                                         genes, ...) {
    indices <- (phyex_set@gene_ids %in% genes)
    
    phyex_set@strata <- phyex_set@strata[indices]
    phyex_set@gene_ids <- phyex_set@gene_ids[indices]
    phyex_set@counts <- phyex_set@counts[indices, ,drop=F]
    
    
    return(phyex_set)
}

#' @title Calculate Transcriptomic Index Confidence Intervals
#' @name TXI_conf_int
#' @description Calculate confidence intervals for TXI values using bootstrap resampling.
#' 
#' @param phyex_set A PhyloExpressionSet object with bootstrapped TXI values
#' @param low_q Lower quantile for confidence interval (default: 0.025)
#' @param high_q Upper quantile for confidence interval (default: 0.975)
#' 
#' @return A list with 'low' and 'high' vectors containing confidence interval bounds
#' 
#' @examples
#' # Calculate 95% confidence intervals
#' # ci <- TXI_conf_int(phyex_set, low_q = 0.025, high_q = 0.975)
#' 
#' @import purrr
#' @importFrom stats quantile
TXI_conf_int <- function(phyex_set, 
                         low_q = .025,
                         high_q = .975) {
    bootstraps <- generate_bootstrapped_txis(phyex_set@pTXI, 
                                             phyex_set@counts_collapsed, 
                                             phyex_set@null_conservation_sample_size)
    CIs <- apply(phyex_set@bootstrapped_txis, 2, stats::quantile, probs=c(low_q, high_q))
    return(list(low=CIs[1, ], high=CIs[2, ]))
}

#' @title Calculate Phylostratum-Specific Transcriptomic Index
#' @description Calculate the phylostratum-specific transcriptomic index (pTXI) for each gene.
#' 
#' @param counts Matrix of expression counts with genes as rows and samples as columns
#' @param strata Numeric or factor vector of phylostratum assignments for each gene
#' 
#' @return Matrix of pTXI values with same dimensions as counts
#' 
#' @details
#' The pTXI is calculated as the relative expression of each gene multiplied by its phylostratum age.
#' Formula: pTXI = (expression / total_expression_per_sample) * phylostratum_age
#' 
#' @examples
#' # Calculate pTXI values
#' # ptxi_values <- pTXI(expression_matrix, phylostratum_vector)
#' 
#' @export
pTXI <- function(counts, strata) {
    
    sweep(counts, 2, colSums(counts), "/") * as.numeric(strata)
}

#' @title Calculate Stratum-Specific Transcriptomic Index
#' @description Calculate the stratum-specific transcriptomic index (sTXI) by summing pTXI values within each phylostratum.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param option Character string specifying calculation method:
#'   - "identity": Sum pTXI values within each stratum
#'   - "add": Cumulative sum across strata
#' 
#' @return Matrix of sTXI values with strata as rows and conditions as columns
#' 
#' @examples
#' # Calculate sTXI values
#' # stxi_values <- sTXI(phyex_set, option = "identity")
#' # stxi_cumsum <- sTXI(phyex_set, option = "add")
#' 
#' @export
sTXI <- function(phyex_set,
                 option = "identity") {
    mat <- rowsum(phyex_set@pTXI, phyex_set@strata)
    if (option == "add") {
        mat <- apply(mat,2,cumsum)
    }
    if (! option %in% c("identity", "add"))
        stop("'option' must be one of 'identity' or 'add",
             call. = FALSE
        )
    return(mat)
}

#' @title Remove Genes from PhyloExpressionSet
#' @description Remove specified genes from a PhyloExpressionSet object.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param genes Character vector of gene IDs to remove
#' @param new_name Character string for the new dataset name (default: auto-generated)
#' 
#' @return A PhyloExpressionSet object with the specified genes removed
#' 
#' @details
#' The returned object shares the same null TXI distribution as the base phyex_set
#' 
#' @examples
#' # Remove specific genes
#' # filtered_set <- remove_genes(phyex_set, c("gene1", "gene2"), 
#' #                             new_name = "Filtered Dataset")
#' 
#' @export
remove_genes <- function(phyex_set, genes, new_name = paste(phyex_set@name, "perturbed")) {
    selected_genes <- setdiff(phyex_set@gene_ids, genes)
    s <- select_genes(phyex_set, selected_genes)
    s@name <- new_name
    s@precomputed_null_conservation_txis <- phyex_set@null_conservation_txis
    return(s)
}

