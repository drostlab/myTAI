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
#' @param strata Factor vector of phylostratum assignments for each gene
#' @param strata_values Numeric vector of phylostratum values used in TXI calculations
#' @param expression Matrix of expression counts with genes as rows and samples as columns
#' @param groups Factor vector indicating which identity each sample belongs to
#' @param name Character string naming the dataset (default: "Phylo Expression Set")
#' @param species Character string specifying the species (default: NULL)
#' @param index_type Character string specifying the transcriptomic index type (default: "TXI")
#' @param identities_label Character string labeling the identities (default: "Identities")
#' @param gene_ids Character vector of gene identifiers (default: character(0), auto-generated from expression rownames if not provided)
#' @param null_conservation_sample_size Numeric value for null conservation sample size (default: 5000)
#' @param .null_conservation_txis Precomputed null conservation TXI values (default: NULL)
#' 
#' @return A PhyloExpressionSetBase object
#' 
#' @details
#' The PhyloExpressionSetBase class serves as the foundation for phylotranscriptomic analysis,
#' providing shared functionality for both bulk and single-cell data types.
#' 
#' \strong{Abstract Properties:}
#' Subclasses must implement the \code{expression_collapsed} property to define how expression
#' data should be collapsed across replicates or cells.
#' 
#' \strong{Computed Properties:}
#' Several properties are computed automatically when accessed:
#' \itemize{
#'   \item \code{gene_ids} - Character vector of gene identifiers (rownames of expression matrix)
#'   \item \code{identities} - Character vector of identity labels (colnames of collapsed expression)
#'   \item \code{sample_names} - Character vector of sample names (colnames of expression matrix)
#'   \item \code{num_identities} - Integer count of unique identities
#'   \item \code{num_samples} - Integer count of total samples
#'   \item \code{num_genes} - Integer count of genes
#'   \item \code{num_strata} - Integer count of phylostrata
#'   \item \code{index_full_name} - Full name of the transcriptomic index type
#'   \item \code{group_map} - List mapping identity names to sample names
#'   \item \code{TXI} - Numeric vector of TXI values for each identity (computed from collapsed expression)
#'   \item \code{TXI_sample} - Numeric vector of TXI values for each sample (computed from raw expression)
#'   \item \code{null_conservation_txis} - Matrix of null conservation TXI values for statistical testing
#' }
#' 
#' \strong{Validation:}
#' The class ensures consistency between expression data, phylostratum assignments, and groupings.
#' All gene-level vectors must have matching lengths, and sample groupings must be consistent.
#' 
#' @import S7
#' @export
PhyloExpressionSetBase <- new_class("PhyloExpressionSetBase",
    properties = list(
        ## REQUIRED PROPERTIES
        strata = new_required_property(
            class = class_factor,
            validator = function(value) {
                if (any(is.na(value))) return("cannot contain NA values. Check phylostratum assignments.")
                if (length(value) == 0) return("cannot be empty. Check phylostratum assignments.")
            },
            name = "strata"
        ),
        strata_values = new_required_property(
            class = class_numeric,
            validator = function(value) {
                if (any(is.na(value))) return("cannot contain NA values. Check phylostratum values.")
                if (length(value) == 0) return("cannot be empty. Check phylostratum values.")
            },
            name = "strata_values"
        ),
        expression = new_required_property(
            # class = class_matrix (S7 doesn't have a class matrix yet)
            validator = function(value) {
                if (any(is.na(value))) return("cannot contain NA values. Check expression matrix.")
                if (length(value) == 0) return("cannot be empty. Check expresison matrix.")
            },
            name = "expression"
        ),  
        groups = new_required_property(
            class = class_factor,
            validator = function(value) {
                if (any(is.na(value))) return("cannot contain NA values. Check groups.")
                if (length(value) == 0) return("cannot be empty. Check groups.")
            },
            name = "groups"
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
        identities_label = new_property(
            class = class_character,
            default = "Identities"
        ),

        ## ABSTRACT PROPERTIES (must be implemented by subclasses)
        expression_collapsed = new_property(
            getter = function(self) stop("expression_collapsed property must be implemented by subclass")
        ),
        
        
        ## SHARED COMPUTED PROPERTIES
        gene_ids = new_property(
            class = class_character,
            getter = function(self) rownames(self@expression),
            setter = function(self, value) {
                if (!length(value))
                    return(self)
                rownames(self@expression) <- value
                self
            },
            validator = function(value) {
                if (any(is.na(value))) return("cannot contain NA values. Check gene IDs.")
                if (length(value) == 0) return("cannot be empty. Check gene IDs.")
                if (length(value) == 1) return("must have at least 2 genes. Check gene IDs.")
            },
            name = "gene_ids"
        ),
        identities = new_property(
            class = class_character,
            getter = function(self) colnames(self@expression_collapsed)
        ),
        sample_names = new_property(
            class = class_character,
            getter = function(self) colnames(self@expression)
        ),
        num_identities = new_property(
            class = class_integer,
            getter = function(self) length(self@identities)
        ),
        num_samples = new_property(
            class = class_integer,
            getter = function(self) length(self@sample_names)
        ),
        num_genes = new_property(
            class = class_integer,
            getter = function(self) length(self@gene_ids)
        ),
        num_strata = new_property(
            class = class_integer,
            getter = function(self) length(levels(self@strata))
        ),
        index_full_name = new_property(
            class = class_character,
            getter = function(self) TI_map[[self@index_type]]
        ),
        group_map = new_property(
            class = class_list,
            getter = function(self) split(self@sample_names, self@groups)
        ),
        
        ## TXI PROPERTIES
        TXI = new_property(
            class = class_double,
            getter = function(self) .TXI(self@expression_collapsed, self@strata_values)
        ),
        TXI_sample = new_property(
            class = class_double,
            getter = function(self) .TXI(self@expression, self@strata_values)
        ),
        
        ## NULL CONSERVATION PROPERTIES
        null_conservation_sample_size = new_property(
            class = class_numeric,
            default = 5000L
        ),
        .null_conservation_txis = new_property(
            default = NULL
        ),
        null_conservation_txis = new_property(
            getter = function(self) {
                if (is.null(self@.null_conservation_txis)) {
                    # Compute and cache the result using expression_collapsed
                    computed_txis <- .memo_generate_conservation_txis(self@strata,
                                                                      self@expression_collapsed,
                                                                      self@null_conservation_sample_size)
                    self@.null_conservation_txis <- computed_txis
                    return(computed_txis)
                } else {
                    return(self@.null_conservation_txis)
                }
            }
        )
    ),
    validator = function(self) {
        if (self@num_genes != length(self@strata)) 
            return("@gene_ids and @strata must have the same length")
        if (length(self@strata) != length(self@strata_values)) 
            return("@strata and @strata_values must have the same length")
        
        # (the above also implicitly checks that the expression matrix and the strata lengths match)
        if (length(self@groups) != self@num_samples) {
            return("@groups must have the same length as the number of columns of @expression")
        }

    }
)




## SHARED FUNCTIONS

#' @title Calculate pTXI for Raw Expression Data
#' @description Internal function to calculate pTXI for expression data.
#' 
#' @param expression_matrix Matrix of expression values
#' @param strata_values Numeric vector of phylostratum values
#' @return Matrix of pTXI values
#' 
#' @keywords internal
.pTXI <- function(expression_matrix, strata_values) {
    relative_expr <- sweep(expression_matrix, 2, colSums(expression_matrix), "/")
    res <- relative_expr * strata_values
    colnames(res) <- colnames(expression_matrix)
    rownames(res) <- rownames(expression_matrix)
    return(res)
}

#' @title Calculate TXI for Raw Expression Data
#' @description Internal function to calculate TXI for expression data.
#' 
#' @param expression_matrix Matrix of expression values
#' @param strata_values Numeric vector of phylostratum values
#' @return Vector of TXI values
#' 
#' @keywords internal
.TXI <- function(expression_matrix, strata_values) {
    return(colSums(.pTXI(expression_matrix, strata_values)))
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

#' @title Collapse PhyloExpressionSet Replicates
#' @description Convert a PhyloExpressionSet with replicates to one with collapsed expression data.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param ... Additional arguments passed to methods
#' 
#' @return A new PhyloExpressionSet object with collapsed expression data
#' 
#' @examples
#' # Collapse replicates in a PhyloExpressionSet
#' collapsed_set <- collapse(example_phyex_set)
#' 
#' @export
collapse <- S7::new_generic("collapse", "phyex_set")

#' @title Select Genes from PhyloExpressionSet
#' @description Extract a subset of genes from a PhyloExpressionSet object.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param ... Additional arguments passed to methods (typically includes 'genes' parameter)
#' 
#' @return A PhyloExpressionSet object containing only the selected genes
#' 
#' @examples
#' # Select specific genes
#' selected_set <- select_genes(example_phyex_set, example_phyex_set@gene_ids[1:10])
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
#' @return Matrix of sTXI values with strata as rows and identities as columns
#' 
#' @examples
#' # Calculate sTXI values
#' stxi_values <- sTXI(example_phyex_set, option = "identity")
#' stxi_cumsum <- sTXI(example_phyex_set, option = "add")
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

#' @title Calculate Phylostratum-Specific Transcriptomic Index
#' @description Calculate pTXI values for expression data. This is a generic function
#' that dispatches based on the input type.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param reps Whether to return pTXI for each sample instead for each group
#' @return Matrix of pTXI values
#' 
#' @export
pTXI <- function(phyex_set, reps=FALSE) {
    if (reps)
        e <- phyex_set@expression
    else
        e <- phyex_set@expression_collapsed
    return(.pTXI(e, phyex_set@strata_values))
}

#' @title Remove Genes from PhyloExpressionSet
#' @description Remove specified genes from a PhyloExpressionSet object.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param genes Character vector of gene IDs to remove
#' @param new_name Character string for the new dataset name (default: auto-generated)
#' @param reuse_null_txis Logical indicating whether to reuse precomputed null conservation TXIs (default: TRUE)
#' 
#' @return A PhyloExpressionSet object with the specified genes removed
#' 
#' @examples
#' # Remove specific genes
#' filtered_set <- remove_genes(example_phyex_set, example_phyex_set@gene_ids[1:5], 
#'                             new_name = "Filtered Dataset")
#' 
#' @export
remove_genes <- function(phyex_set, genes, new_name = paste(phyex_set@name, "perturbed"), reuse_null_txis = TRUE) {
    selected_genes <- setdiff(phyex_set@gene_ids, genes)
    s <- select_genes(phyex_set, selected_genes)
    s@name <- new_name
    
    if (reuse_null_txis)
        s@.null_conservation_txis <- phyex_set@null_conservation_txis
    else
        s@.null_conservation_txis <- NULL
    
    return(s)
}

## BACKWARD COMPATIBILITY FUNCTIONS

#' @title Calculate Transcriptomic Index (TXI)
#' @description Calculate the transcriptomic index values for a PhyloExpressionSet.
#' This function provides backward compatibility with the old TXI() function.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' 
#' @return Numeric vector of TXI values for each identity
#' 
#' @examples
#' # Calculate TXI values
#' txi_values <- TXI(example_phyex_set)
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
#' @return Numeric vector of TAI values for each identity
#' 
#' @examples
#' # Calculate TAI values
#' tai_values <- TAI(example_phyex_set)
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
#' @return Numeric vector of TDI values for each identity
#' 
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
#' @return Numeric vector of TEI values for each identity
#' 
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
#' @return Numeric vector of TPI values for each identity
#' 
#' 
#' @export
TPI <- function(phyex_set) {
    return(phyex_set@TXI)
}

PhyloExpressionSet <- PhyloExpressionSetBase

#' @title Check if object is a PhyloExpressionSet
#' @description Checks if the input is a PhyloExpressionSet S7 object and throws an error if not.
#' @param phyex_set An object to check
#' @return Invisibly returns TRUE if check passes, otherwise throws an error
#' @export
check_PhyloExpressionSet <- function(phyex_set) {
    if (!S7::S7_inherits(phyex_set, PhyloExpressionSetBase)) {
        stop("Input must be a PhyloExpressionSet S7 object.", call. = FALSE)
    }
    invisible(TRUE)
}

#' @title Rename a PhyloExpressionSet
#' @description Returns a copy of the PhyloExpressionSet with a new name.
#' @param phyex_set A PhyloExpressionSet object
#' @param new_name Character string for the new dataset name
#' @return The PhyloExpressionSet object with the updated name
#' @export
rename_phyex_set <- function(phyex_set, new_name) {
    check_PhyloExpressionSet(phyex_set)
    phyex_set@name <- new_name
    phyex_set
}


#' @title Permute Strata in PhyloExpressionSet
#' @description Returns a copy of the PhyloExpressionSet with permuted strata and corresponding strata values.
#' @param phyex_set A PhyloExpressionSet object
#' @return A new PhyloExpressionSet object with permuted strata and strata_values
#' @export
permute_PS <- function(phyex_set) {
    check_PhyloExpressionSet(phyex_set)
    # Permute the strata assignments
    perm_idx <- sample(length(phyex_set@strata))
    new_strata <- phyex_set@strata[perm_idx]
    # Permute the strata_values in the same way
    new_strata_values <- phyex_set@strata_values[perm_idx]
    # Create a copy and assign new strata and strata_values
    new_obj <- phyex_set
    new_obj@strata <- new_strata
    new_obj@strata_values <- new_strata_values
    new_obj
}

#' @title Extract Phylomap from PhyloExpressionSet
#' @description Extract a phylomap (tibble of strata-to-gene mappings) from a PhyloExpressionSet object.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' 
#' @return A tibble with two columns: Stratum (gene age value) and GeneID (gene identifiers)
#' 
#' @examples
#' # Extract phylomap
#' phylomap <- get_phylomap(example_phyex_set)
#' 
#' @export
get_phylomap <- function(phyex_set) {
    check_PhyloExpressionSet(phyex_set)
    tibble::tibble(
        Stratum = phyex_set@strata_values,
        GeneID = phyex_set@gene_ids
    )
}

#' @title Get Strata Legend from PhyloExpressionSet
#' @description Extract the strata legend (phylostratum ranks and names) from a PhyloExpressionSet object.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' 
#' @return A tibble with two columns: Rank (phylostratum values) and Name (phylostratum labels), sorted by Rank
#' 
#' @examples
#' # Get strata legend
#' legend <- get_strata_legend(example_phyex_set)
#' 
#' @export
get_strata_legend <- function(phyex_set) {
    check_PhyloExpressionSet(phyex_set)
    
    tibble::tibble(
        Rank = phyex_set@strata_values,
        Name = as.character(phyex_set@strata)
    ) |>
        dplyr::distinct() |>
        dplyr::arrange(Rank)
}
