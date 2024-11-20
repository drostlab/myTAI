#' @title Test whether GATAI is installed on user machine
#' @description Helper function to test if the user has set up GATAI correctly in their environment.
#' @export
is_installed_gatai <- function() {
  tryCatch({
    gatai <- reticulate::import("gatai", convert = FALSE)
    message("✅ GATAI is installed.")
    return(TRUE)
  }, error = function(e) {
    message("Please check your conda environment and python installation for more details.\n", "Please also consult 'reticulate::install_python()' to install the correct python version to begin with.")
    print(reticulate::py_config())
    stop(simpleError("❌ GATAI is not installed.")) 
  })
}
