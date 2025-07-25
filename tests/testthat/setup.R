data("example_phyex_set")
data("example_phyex_set_old")

# Load single-cell example data if available
example_phyex_set_sc <- NULL
if (requireNamespace("Seurat", quietly = TRUE)) {
    example_phyex_set_sc <- load_example_phyex_set_sc()
}