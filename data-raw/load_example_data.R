
data <- read.csv("data-raw/ath_embryogenesis_2019.csv")
legend <- read.csv("data-raw/strata_legend.tsv", sep="\t")

groups <- rep(c("Preglobular", "Globular", "Early Heart", 
               "Late Heart", "Early Torpedo", "Late Torpedo", 
               "Bent Cotyledon", "Mature Green"), each=3)
example_phyex_set <- as_PhyloExpressionSet(data, 
                                           groups=groups, 
                                           name="Embryogenesis 2019",
                                           species="Arabidopsis thaliana",
                                           index_type="TAI",
                                           strata_labels=legend$Name,
                                           identities_label="Stages")

usethis::use_data(example_phyex_set, overwrite = TRUE)
