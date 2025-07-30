
data <- read.csv("data-raw/ath_embryogenesis_2011.csv", sep="\t")
legend <- read.csv("data-raw/strata_legend.tsv", sep="\t")

example_phyex_set_old <- as_PhyloExpressionSet(data, 
                                               name="Embryogenesis 2011",
                                               species="Arabidopsis thaliana",
                                               index_type="TAI",
                                               strata_legend=legend,
                                               identities_label="Stages")

usethis::use_data(example_phyex_set_old, overwrite = TRUE)
