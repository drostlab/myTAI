
data <- read.csv("data-raw/ath_embryogenesis_2011.csv", sep="\t")
legend <- read.csv("data-raw/strata_legend.tsv", sep="\t")


example_phyex_set_old <- BulkPhyloExpressionSet_from_df(data, 
                                                    name="Embryogenesis 2011",
                                                    species="Arabidopsis thaliana",
                                                    index_type="TAI",
                                                    strata_legend=legend,
                                                    identities_label="Stages")

top_genes <- example_phyex_set_old |> genes_top_mean(p=0.9)
example_phyex_set_old <- example_phyex_set_old |> select_genes(top_genes)

usethis::use_data(example_phyex_set_old, overwrite = TRUE)
