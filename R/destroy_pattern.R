
#' @title Destroy Phylotranscriptomic Pattern Using GATAI
#' @description Apply the GATAI algorithm
#' to identify and remove genes that contribute to phylotranscriptomic patterns.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param num_runs Number of GATAI runs to perform (default: 20)
#' @param runs_threshold Threshold for gene removal consistency across runs (default: 0.5)
#' @param analysis_dir Directory to store GATAI analysis results (default: NULL)
#' @param seed Random seed for reproducibility (default: 1234)
#' 
#' @return A list containing GATAI results including identified genes that contribute to the pattern
#' 
#' @details
#' This function requires the gataiR package to be installed. GATAI systematically removes genes
#' that contribute to phylotranscriptomic patterns by iteratively testing gene removal and
#' evaluating the impact on the overall transcriptomic signature.
#' 
#' @examples
#' # Apply GATAI to identify pattern-contributing genes
#' # gatai_result <- destroy_pattern(phyex_set, num_runs = 10, runs_threshold = 0.6)
#' 
#' @author Filipa Martins Costa
#' @export
destroy_pattern <- function(phyex_set, 
                            num_runs = 20,
                            runs_threshold = 0.5,
                            analysis_dir = NULL,
                            seed = 1234) {
    if (!requireNamespace("gataiR", quietly = TRUE)) {
        stop("Package 'gataiR' must be installed to use this function.")
    }
    
    res <- gataiR::gatai(phyex_set@data_collapsed, 
                         num_runs = num_runs,
                         runs_threshold=runs_threshold, 
                         gataiR_analysis_directory = analysis_dir, 
                         seed = seed)
    
    return(res)
}