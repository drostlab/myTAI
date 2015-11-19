#' @title Perform Phylostratigraphy
#' @description This function takes the protein sequence file of a query organism 
#' and performs the \code{Phylostratigraphy} algorithm to the \code{NCBI nr} database.
#' @param query a character string specifying the path to the protein file of interest (query organism).
#' @param sci.name the scientific name of the query organims.
#' @param seq.type the sequence type of the input file (usually 'protein') for Phylostratigraphy.
#' @param download.databases either a logical in case both \code{NCBI nr} and \code{NCBI Taxonomy} should
#' be loaded automatically to your local machine or a character string specifying which of the both databases
#' shall be downloaded, \code{download.databases = "nr"} (\code{NCBI nr}) or \code{download.databases = "taxdb"} (\code{NCBI Taxonomy})
#' @param max.target.seqs a numeric value specifying the number of aligned sequences to keep.
#' @param eval E-value cutoff for blastp searches. Default is \code{eval = "1E-5"}. 
#' @param cores cores for parallel computations.
#' @param custom.taxa a character vector storing all taxon names for which tax_ids shall be retrieved and gene age assignment based on.
#' @param custom.db a character string specifying the path to the folder in which the customized
#' BLAST database (in FASTA file format) is stored. This customized database is then formatted by \code{makeblastdb} and Phylostratigraphy is performed against this formatted (custom) database.
#' @param nr.path path to the NCBI nr database folder. This can also be specified when \code{download.databases = TRUE} to specify the folder in which the database shall be stored.
#' @param tax.db.path path to the NCBI Taxonomy database folder. This can also be specified when \code{download.databases = TRUE} to specify the folder in which the database shall be stored.
#' @param add.blast.parameters a character string listing the input paramters that shall be passed to the executing BLAST program. Default is \code{NULL}, implicating that a set of default parameters is used when running BLAST.
#' @param add.makeblastdb.params a character string listing the input paramters that shall be passed to the executing makeblastdb program. Default is \code{NULL}, implicating that a set of default parameters is used when running makeblastdb
#' @param blast.path a character string specifying the path to the BLAST program (in case users don't use the default blast path).
#' @author Hajk-Georg Drost
#' @details This
#' @examples 
#' \dontrun{
#' # Construct a phylostratigraphic map for 'Arabidopsis thaliana'
#' Athaliana.PhyloMap <- Phylostratigraphy(
#'    query              = system.file('seqs/example_athal_aa.faa', package = 'orthologr'),
#'    sci.name           = "Arabidopsis thaliana",
#'    seq.type           = "protein", 
#'    download.databases = FALSE, 
#'    max.target.seqs    = 1000000, 
#'    eval               = "1E-5", 
#'    cores              = 1, 
#'    nr.path            = "path/to/local/nr/database/",
#'    tax.db.path        = "path/to/local/Taxonomy/database/")
#' 
#' }
#' @seealso \code{\link{DivergenceStratigraphy}}, \code{\link{TAI}}, \code{\link{RE}}, \code{\link{PlotRE}},
#' \code{\link{PlotPattern}}
#' @import data.table
#' @export

Phylostratigraphy <- function(query,
                              sci.name,
                              seq.type               = "protein",
                              download.databases     = FALSE,
                              max.target.seqs        = 1000000,
                              eval                   = "1E-5",
                              cores                  = 1,
                              custom.taxa            = NULL,
                              custom.db              = NULL,
                              nr.path                = NULL, 
                              tax.db.path            = NULL,
                              add.blast.parameters   = NULL,
                              add.makeblastdb.params = NULL,
                              blast.path             = NULL){
        
        
        if ((!download.databases) & is.null(nr.path) & is.null(tax.db.path))
                stop ("Please specify the path to the location of your local NCBI nr and NCBI Taxonomy databases. Otherwise specify 'download.databases = TRUE' in case these databases shall to be loaded automatically to your local system.")
        
        if ((!is.logical(download.databases)) & (!is.element("nr","taxdb")))
                stop ("Please enter 'download.databases = TRUE' to download the NCBI nr database AND NCBI Taxonomy database or either 'download.databases = 'nr'' for NCBI nr or 'download.databases = 'taxdb'' for NCBI Taxonomy!")
        
        
        # Download NCBI nr and NCBI Taxonomy
        if (download.databases){
                cat("\n")
                
                if (is.null(nr.path))
                        message("The NCBI nr database is loaded to your current working directory: ",getwd()," ...")
                if (!is.null(nr.path))
                        message("The NCBI nr database is loaded to your specified directory: ",nr.path," ...")
                
                cat("\n")
                message("Please make sure that you have constant internet access!")
                cat("\n")
                message("Starting...")
                cat("\n")
                
                if (is.null(nr.path))
                        sapply(biomartr::listDatabases("nr"),biomartr::download_database, path = getwd())
                
                if (!is.null(nr.path))
                        sapply(biomartr::listDatabases("nr"),biomartr::download_database, path = nr.path)
                
                cat("\n")
                
                if (is.null(nr.path))
                        message("The NCBI nr database has successfully been stored in ",getwd()," !")
                
                if (!is.null(nr.path))
                        message("The NCBI nr database has successfully been stored in ",nr.path," !")
                
                cat("\n")
                if (is.null(nr.path)){
                        utils::untar(file.path(getwd(),biomartr::listDatabases("nr")), list = TRUE)
                        utils::untar(file.path(getwd(),biomartr::listDatabases("nr")), exdir = getwd())
                        unlink(file.path(getwd(),biomartr::listDatabases("nr")))
                }
                        
                if (!is.null(nr.path)){
                        utils::untar(file.path(nr.path,biomartr::listDatabases("nr")), list = TRUE)
                        utils::untar(file.path(nr.path,biomartr::listDatabases("nr")), exdir = nr.path)
                        unlink(file.path(nr.path,biomartr::listDatabases("nr")))
                }
                        
                cat("\n")
                
                if (is.null(tax.db.path))
                        message("Now the NCBI Taxonomy database is loaded to your current working directory ",getwd(), "...")
                
                if (!is.null(tax.db.path))
                        message("Now the NCBI Taxonomy database is loaded to your specified directory ",tax.db.path, "...")
                
                cat("\n")
                message("Starting...")
                cat("\n")
                
                if (is.null(tax.db.path)){
                        utils::download.file("ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_prot.zip",file.path(getwd(),"gi_taxid_prot.zip"))
                        utils::download.file("ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip",file.path(getwd(),"taxdmp.zip"))
                        
                        utils::unzip(file.path(getwd(),"gi_taxid_prot.zip"),exdir = getwd())
                        unlink(file.path(getwd(),"gi_taxid_prot.zip"),TRUE,TRUE)
                        utils::unzip(file.path(getwd(),"taxdmp.zip"),exdir = getwd())
                        unlink(file.path(getwd(),"taxdmp.zip"),TRUE,TRUE)
                        
                }
                
                if (!is.null(tax.db.path)){
                        utils::download.file("ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_prot.zip",file.path(tax.db.path,"gi_taxid_prot.zip"))
                        utils::download.file("ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip",file.path(tax.db.path,"taxdmp.zip"))
                        
                        utils::unzip(file.path(tax.db.path,"gi_taxid_prot.zip"),exdir = tax.db.path)
                        unlink(file.path(tax.db.path,"gi_taxid_prot.zip"),TRUE,TRUE)
                        utils::unzip(file.path(tax.db.path,"taxdmp.zip"),exdir = tax.db.path)
                        unlink(file.path(tax.db.path,"taxdmp.zip"),TRUE,TRUE)
                }
                
                cat("\n")
                
                if (is.null(tax.db.path))
                        message("The NCBI Taxonomy database has successfully been stored in ",getwd()," !")
                
                if (!is.null(tax.db.path))
                        message("The NCBI Taxonomy database has successfully been stored in ",tax.db.path," !")
        }
        
        # Download only NCBI nr
        if (download.databases == "nr"){
                cat("\n")
                
                if (is.null(nr.path))
                        message("The NCBI nr database is loaded to your current working directory: ",getwd()," ...")
                if (!is.null(nr.path))
                        message("The NCBI nr database is loaded to your specified directory: ",nr.path," ...")
                
                cat("\n")
                message("Please make sure that you have constant internet access!")
                cat("\n")
                message("Starting...")
                cat("\n")
                
                if (is.null(nr.path))
                        sapply(biomartr::listDatabases("nr"),biomartr::download_database, path = getwd())
                
                if (!is.null(nr.path))
                        sapply(biomartr::listDatabases("nr"),biomartr::download_database, path = nr.path)
                
                cat("\n")
                
                if (is.null(nr.path))
                        message("The NCBI nr database has successfully been stored in ",getwd()," !")
                
                if (!is.null(nr.path))
                        message("The NCBI nr database has successfully been stored in ",nr.path," !")
                
                cat("\n")
        }
        
        # Download only NCBI Taxonomy
        if (download.databases == "taxdb"){
                
                cat("\n")
                
                if (is.null(tax.db.path))
                        message("Now the NCBI Taxonomy database is loaded to your current working directory ",getwd(), "...")
                
                if (!is.null(tax.db.path))
                        message("Now the NCBI Taxonomy database is loaded to your specified directory ",tax.db.path, "...")
                
                cat("\n")
                message("Starting...")
                cat("\n")
                
                if (is.null(tax.db.path)){
                        utils::download.file("ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_prot.zip",file.path(getwd(),"gi_taxid_prot.zip"))
                        utils::download.file("ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip",file.path(getwd(),"taxdmp.zip"))
                }
                
                if (!is.null(tax.db.path)){
                        utils::download.file("ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/gi_taxid_prot.zip",file.path(tax.db.path,"gi_taxid_prot.zip"))
                        utils::download.file("ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip",file.path(tax.db.path,"taxdmp.zip"))
                }
                
                cat("\n")
                
                if (is.null(tax.db.path))
                        message("The NCBI Taxonomy database has successfully been stored in ",getwd()," !")
                
                if (!is.null(tax.db.path))
                        message("The NCBI Taxonomy database has successfully been stored in ",tax.db.path," !")
        }
        
        #query.phylogeny <- taxonomy(organism = sci.name)
        
        if (is.null(custom.db)){
                
                cat("\n")
                cat("\n")
                message("Starting BLAST search against NCBI nr ...")
                cat("\n")
                
                # Perform blastp against the NCBI nr database
                nr.result <- orthologr::blast.nr(query_file      = query,
                                                 nr.path         = nr.path,
                                                 seq_type        = seq.type,
                                                 max.target.seqs = max.target.seqs,
                                                 comp_cores      = cores,
                                                 eval            = eval,
                                                 blast_params    = add.blast.parameters,
                                                 path            = blast.path)
                
                message("Finished BLAST search against NCBI nr !")
                cat("\n")
                
                subject_id <- gi_id <- tax_id <-  parent_tax_id <- sci_name <- NULL
                
                # Parse gi_id's from NCBI nr output 
                nr.result[ , "subject_id"] <-  unlist(lapply(stringr::str_split(nr.result[ , subject_id],"[|]"), function(x) x[2]))
                
                # Read NCBI Taxonomy
                if (is.null(tax.db.path)){
                        cat("\n")
                        message("Reading gi_taxid_prot.dmp file ...")
                        cat("\n")
                        
                        ncbi.taxonomy <- data.table::fread(file.path(getwd(),"gi_taxid_prot.dmp"), sep = "\t")
                        cat("\n")
                        message("Reading nodes.dmp file ...")
                        cat("\n")
                        
                        ncbi.nodes <- readLines(file.path(getwd(),"taxdmp","nodes.dmp"))
                        cat("\n")
                        message("Reading names.dmp file ...")
                        cat("\n")
                        
                        ncbi.names <- readLines(file.path(getwd(),"taxdmp","names.dmp"))
                }
                
                if (!is.null(tax.db.path)){
                        
                        cat("\n")
                        message("Reading gi_taxid_prot.dmp file ...")
                        cat("\n")
                        
                        ncbi.taxonomy <- data.table::fread(file.path(tax.db.path,"gi_taxid_prot.dmp"), sep = "\t")
                        
                        cat("\n")
                        message("Reading nodes.dmp file ...")
                        cat("\n")
                        
                        ncbi.nodes <- readr::read_lines(file.path(tax.db.path,"taxdmp","nodes.dmp"))
                        cat("\n")
                        message("Reading names.dmp file ...")
                        cat("\n")
                        
                        ncbi.names <- readr::read_lines(file.path(tax.db.path,"taxdmp","names.dmp"))
                }
                
                
                # Store NCBI Taxonomy as data.table objects
                data.table::setnames(ncbi.taxonomy, old = c("V1","V2"), new = c("gi_id","tax_id"))
                data.table::setkey(ncbi.taxonomy, gi_id)
                node.dmp <- do.call(rbind,stringr::str_split(ncbi.nodes,"\t[|]\t"))
                colnames(node.dmp) <- c("tax_id","parent_tax_id","rank",
                                        "embl_code","division_id","inherited_div_flag",
                                        "genetic_code_id","inherited_GC_flag", 
                                        "inherited_MGC_flag", "gencode_from_parent", 
                                        "GenBank_hidden_flag", "hidden_subtree root_flag ","comments")
                
                node.dmp.dt <- data.table::data.table(tax_id        = as.numeric(node.dmp[ , 1]), 
                                                      parent_tax_id = as.numeric(node.dmp[ , 2]),
                                                      rank          = node.dmp[ , 3])
                data.table::setkey(node.dmp.dt,tax_id)
                
                names.dmp <- do.call(rbind,stringr::str_split(ncbi.names,"\t[|]\t"))
                colnames(names.dmp) <- c("tax_id","sci_name","uniq_name",
                                        "name_class")
                
                names.dmp.dt <- data.table::data.table(tax_id        = as.numeric(names.dmp[ , 1]), 
                                                       sci_name      = as.character(names.dmp[ , 2]))
                data.table::setkey(names.dmp.dt,sci_name)
                
                cat("\n")
                message("Finished reading!")
                cat("\n")
                
                tax.tree <- vector("numeric")
                tax.tree.sci.names <- vector("character")
                
                cat("\n")
                message("Construct reference taxonomy ...")
                cat("\n")
                
                if (is.null(custom.taxa)){
                        
                        if (is.na(names.dmp.dt[.(sci.name), tax_id]))
                                stop ("Your input reference species '",sci.name,"' could not be found!")
                        # start with the query organism
                        upstream.taxid <- names.dmp.dt[.(sci.name), tax_id]
                        iterator <- 1
                        root <- FALSE
                        
                        while (TRUE){
                                
                                
                                tax.tree[iterator] <- node.dmp.dt[.(upstream.taxid), parent_tax_id]
                                tax.tree.sci.names[iterator] <- names.dmp.dt[tax_id == tax.tree[iterator], sci_name][1]
                                upstream.taxid <- tax.tree[iterator]
                                iterator <- iterator + 1
                                if (upstream.taxid == 131567)
                                        break
                                
                        }
                        
                    Taxonomy.tbl <- data.frame(tax_id = rev(c(names.dmp.dt[.(sci.name), tax_id],tax.tree)),
                                               sci_name = rev(c(names.dmp.dt[.(sci.name), sci_name],tax.tree.sci.names)))
                        
                }
                
                if (!is.null(custom.taxa)){
                        
                        for(i in seq_len(custom.taxa)){
                                
                                tax.tree[i] <- names.dmp.dt[sci_name == custom.taxa[i], tax_id][1]
                        }
                        
                        Taxonomy.tbl <- data.frame(tax_id = rev(tax.tree),
                                                   sci_name = rev(custom.taxa))
                }
                
                cat("\n")
                message("Finished taxonomy construction!")
                cat("\n")
                
                
                
                subj.tbl <- data.table::data.table(subject_id = nr.result[ , subject_id], 
                                       tax_id     = ncbi.taxonomy[.(as.numeric(nr.result[ , subject_id])), tax_id],
                                       parent_tax_id = node.dmp.dt[.(ncbi.taxonomy[.(as.numeric(nr.result[ , subject_id])), tax_id]), parent_tax_id])
                
                nr.res.joined <- dplyr::inner_join(nr.result,subj.tbl, by = "subject_id")
#                 dplyr::group_by(dplyr::group_by(nr.res.joined, query_id),tax_id)
                
                return(nr.res.joined)
                
        }
        
        if (!is.null(custom.db)){
                
                
                if (is.null(add.makeblastdb.params)){
                        system(paste0("makeblastdb -in ",custom.db," -input_type fasta -hash_index ",add.makeblastdb.params))
                }
                
                if (!is.null(add.makeblastdb.params)){
                        system(paste0("makeblastdb -in ",custom.db," -input_type fasta -hash_index ",add.makeblastdb.params))
                }
                
        }
        
        
        
        

        
}




