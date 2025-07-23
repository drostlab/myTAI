data("example_phyex_set")
data("example_phyex_set_old")

# Load single-cell example data if available
if (requireNamespace("Seurat", quietly = TRUE)) {
    data("example_phyex_set_sc")
}

