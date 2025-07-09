
data <- read.csv("data-raw/ath_embryogenesis_2019.csv")
groups <- rep(c("Preglobular", "Globular", "Early Heart", 
               "Late Heart", "Early Torpedo", "Late Torpedo", 
               "Bent Cotyledon", "Mature Green"), each=3)
example_phyex_set <- as_PhyloExpressionSet(data, groups=groups, name="example phyex set")

usethis::use_data(example_phyex_set, overwrite = TRUE)
