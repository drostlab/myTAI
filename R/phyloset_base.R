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

#' @title PhyloExpressionSet Base Class
#' @description Abstract S7 base class for storing and manipulating phylotranscriptomic expression data.
#' This class provides the common interface for both bulk and single-cell phylotranscriptomic data.
#' 
#' @import S7
#' @export
PhyloExpressionSetBase <- new_class("PhyloExpressionSetBase",
    properties = list(
        ## REQUIRED PROPERTIES
        strata = new_required_property(
            class = class_factor,
            validator = function(value) {
                if (any(is.na(value))) "cannot contain NA values. Check phylostratum assignments."
                if (length(value) == 0) "cannot be empty. Check phylostratum assignments."
            },
            name = "strata"
        ),
        gene_ids = new_required_property(
            class = class_character,
            validator = function(value) {
                if (any(is.na(value))) "cannot contain NA values. Check gene IDs."
                if (length(value) == 0) "cannot be empty. Check gene IDs."
            },
            name = "gene_ids"
        ),
        
        ## OPTIONAL PROPERTIES
        name = new_property(
            class = class_character,
            default = "Phylo Expression Set"
        ),
        species = new_property(
            class = class_character,
            default = NULL
        ),
        index_type = new_options_property(
            class = class_character,
            options = names(TI_map),
            default = "TXI"
        ),
        
        ## ABSTRACT PROPERTIES (must be implemented by subclasses)
        identities = new_property(
            class = class_factor,
            getter = function(self) stop("identities property must be implemented by subclass")
        ),
        identities_label = new_property(
            class = class_character,
            default = "Identities"
        ),
        expression_collapsed = new_property(
            getter = function(self) stop("expression_collapsed property must be implemented by subclass")
        ),
        
        ## SHARED COMPUTED PROPERTIES
        num_genes = new_property(
            class = class_integer,
            getter = function(self) as.integer(length(self@gene_ids))
        ),
        num_identities = new_property(
            class = class_integer,
            getter = function(self) as.integer(length(self@identities))
        ),
        num_strata = new_property(
            class = class_integer,
            getter = function(self) as.integer(length(unique(self@strata)))
        ),
        index_full_name = new_property(
            class = class_character,
            getter = function(self) TI_map[[self@index_type]]
        ),
        
        ## TXI PROPERTIES
        TXI = new_property(
            class = class_double,
            getter = function(self) colSums(pTXI(self))
        ),
        TXI_sample = new_property(
            getter = function(self) stop("TXI_sample property must be implemented by subclass")
        ),
        
        ## NULL CONSERVATION PROPERTIES
        null_conservation_sample_size = new_property(
            class = class_numeric,
            default = 5000L
        ),
        precomputed_null_conservation_txis = new_property(
            default = NULL
        ),
        null_conservation_txis = new_property(
            getter = function(self) {
                if (is.null(self@precomputed_null_conservation_txis)) {
                    # Compute and cache the result using expression_collapsed
                    computed_txis <- generate_conservation_txis(self@strata,
                                                              self@expression_collapsed,
                                                              self@null_conservation_sample_size)
                    self@precomputed_null_conservation_txis <- computed_txis
                    return(computed_txis)
                } else {
                    return(self@precomputed_null_conservation_txis)
                }
            }
        ),
        
        ## ABSTRACT SAMPLE NAMES
        sample_names = new_property(
            getter = function(self) stop("sample_names property must be implemented by subclass")
        )
    )
)

## SHARED FUNCTIONS

#' @title Calculate pTXI for Raw Expression Data
#' @description Internal function to calculate pTXI for raw (non-collapsed) expression data.
#' Used for TXI_sample calculations.
#' 
#' @param expression_matrix Matrix of expression values
#' @param strata Vector of phylostratum assignments
#' @return Matrix of pTXI values
#' 
#' @keywords internal
.pTXI_raw <- function(expression_matrix, strata) {
    strata_numeric <- as.numeric(strata)
    relative_expr <- sweep(expression_matrix, 2, colSums(expression_matrix), "/")
    return(relative_expr * strata_numeric)
}

## PRINT METHODS

#' @export
S7::method(print, PhyloExpressionSetBase) <- function(x, ...) {
    cat("PhyloExpressionSet object\n")
    cat("Class:", class(x)[[1]], "\n")
    cat("Name:", x@name, "\n")
    cat("Species:", ifelse(is.null(x@species), "Not specified", x@species), "\n")
    cat("Index type:", x@index_type, "\n")
    cat(x@identities_label, ":", paste(as.character(x@identities), collapse = ", "), "\n")
    cat("Number of genes:", x@num_genes, "\n")
    cat("Number of", tolower(x@identities_label), ":", x@num_identities, "\n")
    cat("Number of phylostrata:", x@num_strata, "\n")
}

## GENERIC FUNCTIONS

#' @title Calculate Phylostratum-Specific Transcriptomic Index
#' @description Calculate pTXI values for expression data. This is a generic function
#' that dispatches based on the input type.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @return Matrix of pTXI values
#' 
#' @export
pTXI <- S7::new_generic("pTXI", dispatch_args = "phyex_set")

#' @export
S7::method(pTXI, PhyloExpressionSetBase) <- function(phyex_set) {
    # Use collapsed expression data for pTXI calculation
    counts <- phyex_set@expression_collapsed
    strata_numeric <- as.numeric(phyex_set@strata)
    
    # Calculate relative expression and multiply by phylostratum
    relative_expr <- sweep(counts, 2, colSums(counts), "/")
    ptxi_matrix <- relative_expr * strata_numeric
    
    # Add proper row and column names
    rownames(ptxi_matrix) <- phyex_set@gene_ids
    colnames(ptxi_matrix) <- as.character(phyex_set@identities)
    
    return(ptxi_matrix)
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
collapse <- S7::new_generic("collapse", "phyex_set")

#' @title Select Genes from PhyloExpressionSet
#' @description Extract a subset of genes from a PhyloExpressionSet object.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param genes Character vector of gene IDs to select
#' 
#' @return A PhyloExpressionSet object containing only the selected genes
#' 
#' @examples
#' # Select specific genes
#' # selected_set <- select_genes(phyex_set, c("gene1", "gene2", "gene3"))
#' 
#' @export
select_genes <- S7::new_generic("select_genes", "phyex_set")

## UTILITY FUNCTIONS

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
    mat <- rowsum(pTXI(phyex_set), phyex_set@strata)
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
#' @examples
#' # Remove specific genes
#' # filtered_set <- remove_genes(phyex_set, c("gene1", "gene2"), 
#' #                             new_name = "Filtered Dataset")
#' 
#' @export
remove_genes <- function(phyex_set, genes, new_name = paste(phyex_set@name, "perturbed"), reuse_null_txis = TRUE) {
    selected_genes <- setdiff(phyex_set@gene_ids, genes)
    s <- select_genes(phyex_set, selected_genes)
    s@name <- new_name
    
    if (reuse_null_txis)
        s@precomputed_null_conservation_txis <- phyex_set@null_conservation_txis
    else
        s@precomputed_null_conservation_txis <- NULL
    
    return(s)
}

## BACKWARD COMPATIBILITY FUNCTIONS

#' @title Calculate Transcriptomic Index (TXI)
#' @description Calculate the transcriptomic index values for a PhyloExpressionSet.
#' This function provides backward compatibility with the old TXI() function.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' 
#' @return Numeric vector of TXI values for each condition
#' 
#' @examples
#' # Calculate TXI values
#' # txi_values <- TXI(phyex_set)
#' 
#' @export
TXI <- function(phyex_set) {
    return(phyex_set@TXI)
}

#' @title Calculate Transcriptomic Age Index (TAI)
#' @description Calculate the transcriptomic age index values for a PhyloExpressionSet.
#' This function provides backward compatibility with the old TAI() function.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' 
#' @return Numeric vector of TAI values for each condition
#' 
#' @examples
#' # Calculate TAI values
#' # tai_values <- TAI(phyex_set)
#' 
#' @export
TAI <- function(phyex_set) {
    return(phyex_set@TXI)
}

#' @title Calculate Transcriptomic Divergence Index (TDI)
#' @description Calculate the transcriptomic divergence index values for a PhyloExpressionSet.
#' This function provides backward compatibility with the old TDI() function.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' 
#' @return Numeric vector of TDI values for each condition
#' 
#' @examples
#' # Calculate TDI values
#' # tdi_values <- TDI(phyex_set)
#' 
#' @export
TDI <- function(phyex_set) {
    return(phyex_set@TXI)
}

#' @title Calculate Transcriptomic Evolutionary Index (TEI)
#' @description Calculate the transcriptomic evolutionary index values for a PhyloExpressionSet.
#' This function provides backward compatibility with the old TEI() function.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' 
#' @return Numeric vector of TEI values for each condition
#' 
#' @examples
#' # Calculate TEI values
#' # tei_values <- TEI(phyex_set)
#' 
#' @export
TEI <- function(phyex_set) {
    return(phyex_set@TXI)
}

#' @title Calculate Transcriptomic Polymorphism Index (TPI)
#' @description Calculate the transcriptomic polymorphism index values for a PhyloExpressionSet.
#' This function provides backward compatibility with the old TPI() function.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' 
#' @return Numeric vector of TPI values for each condition
#' 
#' @examples
#' # Calculate TPI values
#' # tpi_values <- TPI(phyex_set)
#' 
#' @export
TPI <- function(phyex_set) {
    return(phyex_set@TXI)
}
