data("example_phyex_set")
data("example_phyex_set_old")

# Load single-cell data if available
# This will be added later by the user
try({
    data("example_phyex_set_sc")
}, silent = TRUE)

set.seed(1234)