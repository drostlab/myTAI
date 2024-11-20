#' @title Low level interface function with GATAI
#' @description This function interfaces with the GATAI software which removes the genes that create the \code{\link{PlotSignature}} pattern.
#' @param ExpressionSet a standard PhyloExpressionSet, DivergenceExpressionSet or PolymorphismsExpressionSet object.
#' @param singlecell_data a logical value indicating whether or not the input \code{ExpressionSet} is a single cell dataset.
#' @param gatai_results_path where shall GATAI results be stored. Default \code{gatai_results_path = tempdir()}.
#' @param result_overwrite Logical value indicating whether local results of previous GATAI run should be over-written and re-computed?
#' @param condaenv The conda environment to use. For \code{\link[reticulate]{use_condaenv}}, this can be the name, the absolute prefix path, or the absolute path to the python binary. If the name is ambiguous, the first environment is used and a warning is issued. For \code{\link[reticulate]{use_miniconda}}, the only conda installation searched is the one installed by \code{\link[reticulate]{install_miniconda}}.
#' @author Filipa Martins Costa and Hajk-Georg Drost
#' @export
GATAI <- function(ExpressionSet,
                  singlecell_data = FALSE,
                  gatai_results_path = tempdir(),
                  result_overwrite = FALSE,
                  condaenv = NULL) {
  
  # use specific conda environment if selected
  if (!is.null(condaenv)) {
    reticulate::use_condaenv(condaenv = condaenv, required = TRUE)
  }
  
  # check if GATAI is properly installed
  is_installed_gatai()
  message("Starting GATAI optimisation ...")
  
  # generate a tempdir version of the input ExpressionSet to run GATAI with random id attached
  export_expressionset_for_gatai <- file.path(tempdir(), "expressionset.tsv")

  # file path to extracted genes output file
  extracted_genes_path <- file.path(gatai_results_path, "removed_genes", "extracted_genes.txt")
  
  # check if a local version of GATAI result is already present
  if(fs::file_exists(extracted_genes_path) && !result_overwrite){
    message("It seems like results from a previous GATAI run were still present at '", extracted_genes_path, "'.")
    message("This previous GATAI result was now imported. Please select argument 'result_overwrite = TRUE' in case you would like to re-calculate eveything with GATAI.")
    removed_genes <- readr::read_tsv(
      file.path(gatai_results_path, "removed_genes", "extracted_genes.txt"),
      col_names = F,
      show_col_types = FALSE
    )
    
    colnames(removed_genes) <- "removed_genes"
  } else {
    message(
      "A local version of the input ExpressionSet was stored for GATAI runs at '",
      export_expressionset_for_gatai,
      "'."
    )
    
    # store input ExpressionSet in tempdir()
    readr::write_delim(
      ExpressionSet,
      file = export_expressionset_for_gatai,
      delim = "\t",
      col_names = TRUE
    )
    
    # run GATAI
    gatai <- reticulate::import("gatai")
    
    if (!singlecell_data) {
      # here we need to specify a path to tempdir()/gatai_results_path where removed genes are stored
      gatai_command <- paste(
        "gatai run_minimizer",
        export_expressionset_for_gatai,
        file.path(gatai_results_path, "removed_genes"),
        "--save_stats"
      )
      
    } else {
      # here we need to specify a path to tempdir()/gatai_results_path where removed genes are stored
      gatai_command <- paste(
        "gatai run_minimizer",
        export_expressionset_for_gatai,
        file.path(gatai_results_path,
                  "removed_genes")
        , "--single_cell", "--save_stats")
    }
    
    # systems call with tryCatch
    tryCatch(
          {
            system(gatai_command)
            message("GATAI Optimization Complete.")
          },
          warning = function(w) {
            message("Warning during GATAI optimization: ", conditionMessage(w))
          },
          error = function(e) {
            message("Error during GATAI optimization: ", conditionMessage(e))
          }
        )
    
    # print the path to GATAI results as message  
    message(
      "GATAI results are stored at '",
      file.path(gatai_results_path, "removed_genes"),
      "'."
    )
    
    # read the extracted_genes.txt file
    removed_genes <- readr::read_tsv(
      file.path(gatai_results_path, "removed_genes", "extracted_genes.txt"),
      col_names = F,
      show_col_types = FALSE
    )
    
    colnames(removed_genes) <- "removed_genes"
    
  }

  # read the summary statistics file
  statistics <- readr::read_tsv(
      file.path(gatai_results_path, "removed_genes", "summary.txt"),
      col_names = F,
      show_col_types = FALSE
  )
    
  # print summary statistics as messages
  colnames(statistics) <- "statistics"
  message("Summary statistics:")
  for (summary in statistics$statistics) {
      message(summary)
    }

  
  return(removed_genes)
}
