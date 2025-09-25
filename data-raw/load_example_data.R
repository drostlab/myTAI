library(myTAI)

data <- read.csv("data-raw/ath_embryogenesis_2019.csv")
legend <- read.csv("data-raw/strata_legend.tsv", sep="\t")

groups <- rep(c("Preglobular", "Globular", "Early Heart", 
               "Late Heart", "Early Torpedo", "Late Torpedo", 
               "Bent Cotyledon", "Mature Green"), each=3)
example_phyex_set <- BulkPhyloExpressionSet_from_df(data, 
                                                    groups=groups, 
                                                    name="Embryogenesis 2019",
                                                    species="Arabidopsis thaliana",
                                                    index_type="TAI",
                                                    strata_legend=legend,
                                                    identities_label="Stages") |>
                                                    remove_genes(c("AT2G01021"), 
                                                    new_name="Embryogenesis 2019",
                                                    reuse_null_txis=FALSE)

top_genes <- example_phyex_set |> genes_top_mean(p=0.93)
example_phyex_set <- example_phyex_set |> select_genes(top_genes)                                                

usethis::use_data(example_phyex_set, overwrite = TRUE)
