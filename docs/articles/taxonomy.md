# Obtaining taxonomy information

The `taxonomy()` function (formerly) implemented in `myTAI` relies on
the powerful package [taxize](https://github.com/ropensci/taxize). More
specifically, the taxonomic information retrieval has been customized
for the `myTAI` standard and for organism specific information
retrieval.

While the previous `taxonomy()` function has been **deprecated** since
`taxize` was pulled from CRAN, users can nevertheless follow the
taxonomy pipeline by installing the `taxize` package and **copy** the
old taxonomy function.

``` r
# install taxize from CRAN
install.packages("taxize")

# if taxize is not available again
install.packages("remotes")
remotes::install_github("ropensci/taxize")
```

**Copy** the taxonomy function:

**open for the taxonomy function**

Click on the copy icon to copy the function.

``` r
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
```

The `taxonomy()` function can be used to classify genomes according to
phylogenetic classification into Phylostrata (Phylostratigraphy) or to
retrieve species specific taxonomic information when performing
Divergence Stratigraphy.

For larger taxonomy queries it may be useful to create an NCBI Account
and set up an [ENTREZ API
KEY](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/).

``` r
# install.packages(c("taxize", "usethis"))
taxize::use_entrez()
# Create your key from your (brand-new) account's. 
# After generating your key set it as ENTREZ_KEY in .Renviron.
# ENTREZ_KEY='youractualkeynotthisstring'
# For that, use usethis::edit_r_environ()
usethis::edit_r_environ()
```

## Taxonomic Information Retrieval

The `taxonomy()` function to retrieve taxonomic information.

**retrieve taxonomy hierarchy**

In the following example we will obtain the taxonomic hierarchy of
`Arabidopsis thaliana` from [NCBI
Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy).

``` r
# retrieving the taxonomic hierarchy of "Arabidopsis thaliana"
# from NCBI Taxonomy
taxonomy( organism = "Arabidopsis thaliana", 
          db       = "ncbi",
          output   = "classification" )
```

Show output

    ## ══  1 queries  ═══════════════

    ## ✔  Found:  Arabidopsis+thaliana
    ## ══  Results  ═════════════════
    ## 
    ## • Total: 1 
    ## • Found: 1 
    ## • Not Found: 0

    ##                    name          rank      id
    ## 1    cellular organisms cellular root  131567
    ## 2             Eukaryota        domain    2759
    ## 3         Viridiplantae       kingdom   33090
    ## 4          Streptophyta        phylum   35493
    ## 5        Streptophytina     subphylum  131221
    ## 6           Embryophyta         clade    3193
    ## 7          Tracheophyta         clade   58023
    ## 8         Euphyllophyta         clade   78536
    ## 9         Spermatophyta         clade   58024
    ## 10        Magnoliopsida         class    3398
    ## 11      Mesangiospermae         clade 1437183
    ## 12       eudicotyledons         clade   71240
    ## 13           Gunneridae         clade   91827
    ## 14         Pentapetalae         clade 1437201
    ## 15               rosids         clade   71275
    ## 16              malvids         clade   91836
    ## 17          Brassicales         order    3699
    ## 18         Brassicaceae        family    3700
    ## 19           Camelineae         tribe  980083
    ## 20          Arabidopsis         genus    3701
    ## 21 Arabidopsis thaliana       species    3702

The `organism` argument takes the scientific name of a query organism,
the `db` argument specifies that database from which the corresponding
taxonomic information shall be retrieved, e.g. `ncbi` (NCBI Taxonomy)
and `itis` (Integrated Taxonomic Information System) and the `output`
argument specifies the type of taxonomic information that shall be
returned for the query organism, e.g. `classification`, `taxid`, or
`children`.

The output of `classification` is a `data.frame` storing the taxonomic
hierarchy of `Arabidopsis thaliana` starting with `cellular organisms`
up to `Arabidopsis thaliana`. The first column stores the taxonomic
name, the second column the taxonomic rank, and the third column the
NCBI Taxonomy id for corresponding taxa.

Analogous `classification` information can be obtained from different
databases.

``` r
# retrieving the taxonomic hierarchy of "Arabidopsis thaliana"
# from the Integrated Taxonomic Information System
taxonomy( organism = "Arabidopsis thaliana", 
          db       = "itis",
          output   = "classification" )
```

Show output

    ## ══  1 queries  ═══════════════

    ## ✔  Found:  Arabidopsis thaliana
    ## ══  Results  ═════════════════
    ## 
    ## • Total: 1 
    ## • Found: 1 
    ## • Not Found: 0

    ##                    name          rank     id
    ## 1               Plantae       kingdom 202422
    ## 2         Viridiplantae    subkingdom 954898
    ## 3          Streptophyta  infrakingdom 846494
    ## 4           Embryophyta superdivision 954900
    ## 5          Tracheophyta      division 846496
    ## 6       Spermatophytina   subdivision 846504
    ## 7         Magnoliopsida         class  18063
    ## 8               Rosanae    superorder 846548
    ## 9           Brassicales         order 822943
    ## 10         Brassicaceae        family  22669
    ## 11          Arabidopsis         genus  23040
    ## 12 Arabidopsis thaliana       species  23041

The `output` argument allows you to directly access taxonomy ids for a
query organism or species.

**retrieve taxonomy ID from `ncbi`**

``` r
# retrieving the taxonomy id of the query organism from NCBI Taxonomy
taxonomy( organism = "Arabidopsis thaliana", 
          db       = "ncbi", 
          output   = "taxid" )
```

Show output

    ## ══  1 queries  ═══════════════

    ## ✔  Found:  Arabidopsis+thaliana
    ## ══  Results  ═════════════════
    ## 
    ## • Total: 1 
    ## • Found: 1 
    ## • Not Found: 0

    ##     id
    ## 1 3702

**retrieve taxonomy ID from `itis`**

``` r
# retrieving the taxonomy id of the query organism from Integrated Taxonomic Information Service
taxonomy( organism = "Arabidopsis", 
          db       = "itis", 
          output   = "taxid" )
```

Show output

    ## ══  1 queries  ═══════════════

    ## ✔  Found:  Arabidopsis
    ## ══  Results  ═════════════════
    ## 
    ## • Total: 1 
    ## • Found: 1 
    ## • Not Found: 0

    ##      id
    ## 1 23040

So far, the following data bases can be accesses to retrieve taxonomic
information:

- `db = "itis"` : Integrated Taxonomic Information Service
- `db = "ncbi"` : National Center for Biotechnology Information

**How does the `taxonomy(db = "ncbi")` output differ from `GenEra`?**

The taxonomic classifications should be the same between
`taxonomy(..., db = "ncbi")` and the taxonomic classifications in the
`GenEra` output (since it uses NCBI taxdump as input). But it should be
noted that the recent updates to NCBI taxonomy has meant that the
highest order ranks (cellular root, domain, kingdom etc.) may differ.

## Retrieve Children Nodes

Another `output` supported by `taxonomy()` is `children` that returns
the immediate children taxa for a query organism. This feature is useful
to determine species relationships for quantifying recent evolutionary
conservation with Divergence Stratigraphy.

**retrieve children nodes from `ncbi`**

``` r
# retrieve children taxa of the query organism stored in the correspondning database
taxonomy( organism = "Arabidopsis", 
          db       = "ncbi", 
          output   = "children" )
```

Show output

    ## ══  1 queries  ═══════════════

    ## ✔  Found:  Arabidopsis
    ## ══  Results  ═════════════════
    ## 
    ## • Total: 1 
    ## • Found: 1 
    ## • Not Found: 0

    ##    childtaxa_id
    ## 1       2608267
    ## 2       2486701
    ## 3       1837063
    ## 4       1547872
    ## 5       1328956
    ## 6       1240361
    ## 7        869751
    ## 8        869750
    ## 9        864766
    ## 10       412662
    ## 11       378006
    ## 12       347883
    ## 13       302551
    ## 14        97980
    ## 15        97979
    ## 16        81970
    ## 17        59690
    ## 18        59689
    ## 19        45251
    ## 20        45249
    ## 21        38785
    ## 22         3702
    ##                                                        childtaxa_name
    ## 1                            Arabidopsis arenosa x Arabidopsis lyrata
    ## 2                            Arabidopsis lyrata x Arabidopsis halleri
    ## 3                          Arabidopsis thaliana x Arabidopsis halleri
    ## 4                                               Arabidopsis umezawana
    ## 5  (Arabidopsis thaliana x Arabidopsis arenosa) x Arabidopsis suecica
    ## 6                          Arabidopsis thaliana x Arabidopsis arenosa
    ## 7         Arabidopsis thaliana x Arabidopsis halleri subsp. gemmifera
    ## 8                           Arabidopsis thaliana x Arabidopsis lyrata
    ## 9                                         Arabidopsis septentrionalis
    ## 10                                            Arabidopsis pedemontana
    ## 11                         Arabidopsis arenosa x Arabidopsis thaliana
    ## 12                                              Arabidopsis arenicola
    ## 13                                              Arabidopsis petrogena
    ## 14                                               Arabidopsis croatica
    ## 15                                            Arabidopsis cebennensis
    ## 16                                                Arabidopsis halleri
    ## 17                                             Arabidopsis kamchatica
    ## 18                                                 Arabidopsis lyrata
    ## 19                                               Arabidopsis neglecta
    ## 20                                                Arabidopsis suecica
    ## 21                                                Arabidopsis arenosa
    ## 22                                               Arabidopsis thaliana
    ##    childtaxa_rank
    ## 1         species
    ## 2         species
    ## 3         species
    ## 4         species
    ## 5         species
    ## 6         species
    ## 7         species
    ## 8         species
    ## 9         species
    ## 10        species
    ## 11        species
    ## 12        species
    ## 13        species
    ## 14        species
    ## 15        species
    ## 16        species
    ## 17        species
    ## 18        species
    ## 19        species
    ## 20        species
    ## 21        species
    ## 22        species

**retrieve children nodes from `itis`**

``` r
# retrieve children taxa of the query organism stored in the correspondning database
taxonomy( organism = "Arabidopsis", 
          db       = "itis", 
          output   = "children" )
```

Show output

    ## ══  1 queries  ═══════════════

    ## ✔  Found:  Arabidopsis
    ## ══  Results  ═════════════════
    ## 
    ## • Total: 1 
    ## • Found: 1 
    ## • Not Found: 0
    ## ══  1 queries  ═══════════════

    ## ✔  Found:  Arabidopsis
    ## ══  Results  ═════════════════
    ## 
    ## • Total: 1 
    ## • Found: 1 
    ## • Not Found: 0

    ##    parentname parenttsn rankname             taxonname    tsn
    ## 1 Arabidopsis     23040  Species  Arabidopsis thaliana  23041
    ## 2 Arabidopsis     23040  Species   Arabidopsis arenosa 823130
    ## 3 Arabidopsis     23040  Species    Arabidopsis lyrata 823171
    ## 4 Arabidopsis     23040  Species Arabidopsis arenicola 823113

These results allow us to choose `subject` organisms for [Divergence
Stratigraphy](https://drostlab.github.io/myTAI/articles/other-strata.md).
