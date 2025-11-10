#' @title Retrieve taxonomy categories from NCBI Taxonomy
#' @description This function retrieves category information 
#' from NCBI Taxonomy and is able to filter kingdom specific taxids.
#' @param db.path path to download and store the NCBI Taxonomy \code{categories.dmp} file. Default is the \code{tempdir()}
#' directory.
#' @param download a logical value specifying whether or not the \code{categories.dmp} shall be downloaded (\code{download = TRUE})
#' or whether a local version already exists on the users machine (\code{download = TRUE} - in this case please specify the \code{db.path} argument to target the local \code{categories.dmp} file).
#' @param update should the local file be updated? Please specify the \code{db.path} argument to target the local \code{categories.dmp} file.
#' @param filter a character string specifying the kingdom of life for which taxids shall be returned. Options are
#' \code{"Archea"}, \code{"Bacteria"}, \code{"Eukaryota"}, \code{"Viruses"}, \code{"Unclassified"}.
#' @author Hajk-Georg Drost
#' 
#' @return A tibble object containing taxonomy category information
#'
#' @examplesIf FALSE 
#' # download categories.dmp file to current working directory 
#' # and filter for 'Archea' taxids
#' Archea.taxids <- taxid(db.path = getwd(), filter = "Archea", download = TRUE)
#' 
#' # Once the NCBI Taxonomy 'categories.dmp' file is downloaded to your machine ('download = TRUE')
#' # the 'taxid()' function can be proceed on the local 'categories.dmp' file
#' # e.g. filter for Virus taxids
#' Virus.taxids <- taxid(db.path = getwd(), filter = "Viruses")
#' 
#' @export


taxid <- function(db.path, download = FALSE, update = FALSE, filter = NULL){
        
        if (download || update){
                tryCatch({
                        utils::download.file("ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxcat.zip",file.path(db.path,"taxcat.zip"))
                        utils::unzip(file.path(db.path,"taxcat.zip"))
                }, error = function(e) {
                        warning("Failed to download taxonomy data. Please check your internet connection or try again later.", call. = FALSE)
                        return(NULL)
                })
        }
        
        top_level_category <- NULL
        
        if (!file.exists(file.path(db.path, "categories.dmp"))) {
                warning("categories.dmp file not found. Please download it first using download = TRUE.", call. = FALSE)
                return(NULL)
        }
        
        ncbi.taxcat <- tryCatch({
                readr::read_tsv(file.path(db.path, "categories.dmp"), col_names = TRUE, show_col_types = FALSE)
        }, error = function(e) {
                warning("Failed to read categories.dmp file. The file may be corrupted.", call. = FALSE)
                return(NULL)
        })
        
        if (is.null(ncbi.taxcat)) {
                return(NULL)
        }
        
        colnames(ncbi.taxcat) <- c("top_level_category","species_level_taxid","tax_id")

        if (is.null(filter))
                return (ncbi.taxcat)
        
        if (!is.null(filter)){
                
                if (!is.element(filter, c("Archea","Bacteria","Eukaryota","Viruses","Unclassified")))
                        stop ("Please choose from the following filters: Archea, Bacteria , Eukaryota , Viruses , Unclassified")
                
                if (is.element(filter, c("Archea","Bacteria","Eukaryota","Viruses","Unclassified"))){
                        
                     filtered.taxids <-  switch(filter, 
                                                Archea       = dplyr::filter(dplyr::tbl_df(ncbi.taxcat),top_level_category == "A"),
                                                Bacteria     = dplyr::filter(dplyr::tbl_df(ncbi.taxcat),top_level_category == "B"),
                                                Eukaryota    = dplyr::filter(dplyr::tbl_df(ncbi.taxcat),top_level_category == "E"),
                                                Viruses      = dplyr::filter(dplyr::tbl_df(ncbi.taxcat),top_level_category == "V"),
                                                Unclassified = dplyr::filter(dplyr::tbl_df(ncbi.taxcat),top_level_category == "U")
                               )
                }
                
                if (nrow(filtered.taxids) == 0) {
                        warning("No taxids found for filter '", filter, "'.", call. = FALSE)
                        return(NULL)
                }
                
                return(filtered.taxids)
        }
}





