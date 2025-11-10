## ----message=FALSE, warning=FALSE, results='hide'-----------------------------
#' @title Retrieving Taxonomic Information of a Query Organism
#' @description This function takes the scientific name of a query organism
#' and returns selected output formats of taxonomic information for the corresponding organism.
#' @param organism a character string specifying the scientific name of a query organism.
#' @param db a character string specifying the database to query, e.g. \code{db} = \code{"itis"} or \code{"ncbi"}.
#' @param output a character string specifying the taxonomic information that shall be returned. 
#' Implemented are: \code{output} = \code{"classification"}, \code{"taxid"}, or \code{"children"}.
#' @details This function is based on the powerful package \pkg{taxize} and implements
#' the customized retrieval of taxonomic information for a query organism. 
#' 
#' The following data bases can be selected to retrieve taxonomic information:
#' 
#' \itemize{
#' \item \code{db = "itis"} : Integrated Taxonomic Information Service
#' \item \code{db = "ncbi"} : National Center for Biotechnology Information
#' }
#' 
#' 
#' 
#' @author Hajk-Georg Drost
#' @examples
#' \dontrun{
#' # retrieving the taxonomic hierarchy of "Arabidopsis thaliana"
#' # from NCBI Taxonomy
#' taxonomy("Arabidopsis thaliana",db = "ncbi")
#' 
#' # the same can be applied to database : "itis"
#'  taxonomy("Arabidopsis thaliana",db = "itis")
#' 
#' # retrieving the taxonomic hierarchy of "Arabidopsis"
#'  taxonomy("Arabidopsis",db = "ncbi") # analogous : db = "ncbi" or "itis"
#' 
#' # or just "Arabidopsis"
#'  taxonomy("Arabidopsis",db = "ncbi")
#' 
#' # retrieving the taxonomy id of the query organism and in the correspondning database
#' # taxonomy("Arabidopsis thaliana",db = "ncbi", output = "taxid")
#' 
#' # the same can be applied to databases : "ncbi" and "itis"
#'  taxonomy("Arabidopsis thaliana",db = "ncbi", output = "taxid")
#'  taxonomy("Arabidopsis thaliana",db = "itis", output = "taxid")
#' 
#' 
#' # retrieve children taxa of the query organism stored in the correspondning database
#'  taxonomy("Arabidopsis",db = "ncbi", output = "children")
#' 
#' # the same can be applied to databases : "ncbi" and "itis"
#'  taxonomy("Arabidopsis thaliana",db = "ncbi", output = "children")
#'  taxonomy("Arabidopsis thaliana",db = "itis", output = "children")
#'  
#' }
#' @references
#' 
#' Scott Chamberlain and Eduard Szocs (2013). taxize - taxonomic search and retrieval in R. F1000Research,
#' 2:191. URL: http://f1000research.com/articles/2-191/v2.
#' 
#' Scott Chamberlain, Eduard Szocs, Carl Boettiger, Karthik Ram, Ignasi Bartomeus, and John Baumgartner
#' (2014) taxize: Taxonomic information from around the web. R package version 0.3.0.
#' https://github.com/ropensci/taxize
#' @export

taxonomy <- function(organism, db = "ncbi", output = "classification"){
        
        if (!is.element(output,c("classification","taxid","children")))
                stop ("The output '",output,"' is not supported by this function.")
        
        if (!is.element(db,c("ncbi","itis")))
                stop ("Database '",db,"' is not supported by this function.")
        
        name <- id <- NULL

        tax_hierarchy <- tryCatch({
                if (db == "ncbi")
                        as.data.frame(taxize::classification(taxize::get_uid(organism), db = "ncbi")[[1]])
                else if (db == "itis")    
                        as.data.frame(taxize::classification(taxize::get_tsn(organism), db = "itis")[[1]])
        }, error = function(e) {
                warning("Could not retrieve taxonomy information from ", db, ". Check internet connection or try again later.", call. = FALSE)
                return(NULL)
        })
        
        if (is.null(tax_hierarchy)) {
                return(NULL)
        }
        
        if(output == "classification"){
                
                return(tax_hierarchy)
        }
        
        if(output == "taxid"){
                
                        return(dplyr::select(dplyr::filter(tax_hierarchy, name == organism),id))
        }
        
        if(output == "children"){
                
                return(as.data.frame(taxize::children(organism, db = db)[[1]]))
        } 
}

## ----eval=FALSE---------------------------------------------------------------
# # retrieving the taxonomic hierarchy of "Arabidopsis thaliana"
# # from NCBI Taxonomy
# taxonomy( organism = "Arabidopsis thaliana",
#           db       = "ncbi",
#           output   = "classification" )

## ----message = FALSE, warning = FALSE, echo = FALSE, eval = requireNamespace("taxize", quietly = TRUE) && curl::has_internet()----
taxonomy( organism = "Arabidopsis thaliana", 
          db       = "ncbi",
          output   = "classification" )

## ----eval=FALSE---------------------------------------------------------------
# # retrieving the taxonomic hierarchy of "Arabidopsis thaliana"
# # from the Integrated Taxonomic Information System
# taxonomy( organism = "Arabidopsis thaliana",
#           db       = "itis",
#           output   = "classification" )

## ----message = FALSE, warning = FALSE, echo = FALSE, eval = requireNamespace("taxize", quietly = TRUE) && curl::has_internet()----
taxonomy( organism = "Arabidopsis thaliana", 
          db       = "itis",
          output   = "classification" )

## ----eval=FALSE---------------------------------------------------------------
# # retrieving the taxonomy id of the query organism from NCBI Taxonomy
# taxonomy( organism = "Arabidopsis thaliana",
#           db       = "ncbi",
#           output   = "taxid" )

## ----message = FALSE, warning = FALSE, echo = FALSE, eval = requireNamespace("taxize", quietly = TRUE) && curl::has_internet()----
taxonomy( organism = "Arabidopsis thaliana", 
          db       = "ncbi", 
          output   = "taxid" )

## ----eval=FALSE---------------------------------------------------------------
# # retrieving the taxonomy id of the query organism from Integrated Taxonomic Information Service
# taxonomy( organism = "Arabidopsis",
#           db       = "itis",
#           output   = "taxid" )

## ----message = FALSE, warning = FALSE, echo = FALSE, eval = requireNamespace("taxize", quietly = TRUE) && curl::has_internet()----
taxonomy( organism = "Arabidopsis", 
          db       = "itis", 
          output   = "taxid" )

## ----eval=FALSE---------------------------------------------------------------
# # retrieve children taxa of the query organism stored in the correspondning database
# taxonomy( organism = "Arabidopsis",
#           db       = "ncbi",
#           output   = "children" )

## ----message = FALSE, warning = FALSE, echo = FALSE, eval = requireNamespace("taxize", quietly = TRUE) && curl::has_internet()----
taxonomy( organism = "Arabidopsis", 
          db       = "ncbi", 
          output   = "children" )

## ----eval=FALSE---------------------------------------------------------------
# # retrieve children taxa of the query organism stored in the correspondning database
# taxonomy( organism = "Arabidopsis",
#           db       = "itis",
#           output   = "children" )

## ----message = FALSE, warning = FALSE, echo = FALSE, eval = requireNamespace("taxize", quietly = TRUE) && curl::has_internet()----
taxonomy( organism = "Arabidopsis", 
          db       = "itis", 
          output   = "children" )

