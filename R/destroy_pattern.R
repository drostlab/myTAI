
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