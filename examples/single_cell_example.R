# Single-Cell PhyloExpressionSet Example
# This example demonstrates how to use the ScPhyloExpressionSet class
# for single-cell phylotranscriptomic analysis

# Required packages
library(myTAI)  # This package
library(Seurat) # For single-cell data handling
library(ggplot2) # For plotting

# Example workflow:

# 1. Load your Seurat object and phylostratum mapping
# seurat_obj <- readRDS("path/to/your/seurat_object.rds")
# phylo_map <- read.csv("path/to/phylostratum_map.csv")

# 2. Create a ScPhyloExpressionSet
# Assumes your Seurat object has cell type annotations in "celltype" column
# sc_phyex_set <- as_ScPhyloExpressionSet(
#     seurat = seurat_obj,
#     phylomap = phylo_map,
#     cell_identity = "celltype",
#     slot = "data",
#     name = "Brain Development Dataset",
#     min_cells_per_type = 20
# )

# 3. Explore the data
# print(sc_phyex_set)
# summary(sc_phyex_set@cell_tai)

# 4. Plotting functions

# Violin plot showing TAI distribution across cell types
# tai_violin <- plot_sc_signature(sc_phyex_set, violin = TRUE, points = FALSE)
# print(tai_violin)

# UMAP plot colored by cell type
# umap_celltype <- plot_sc_sample_space(sc_phyex_set, "umap", "celltype")
# print(umap_celltype)

# UMAP plot colored by TAI values
# umap_tai <- plot_sc_sample_space(sc_phyex_set, "umap", "TAI")
# print(umap_tai)

# Summary plot with error bars
# tai_summary <- plot_sc_tai_summary(sc_phyex_set, error_type = "se")
# print(tai_summary)

# 5. Statistical analysis

# Compare TAI across cell types
# comparison <- compare_sc_tai(sc_phyex_set, test = "anova", pairwise = TRUE)
# print(comparison$overall_p)
# print(comparison$pairwise)

# 6. Data extraction

# Get cell-level data with coordinates
# cell_data <- get_cell_tai_data(sc_phyex_set, include_coords = TRUE, reduction = "umap")
# head(cell_data)

# 7. Subsetting

# Subset to specific cell types
# subset_set <- subset_cell_types(sc_phyex_set, c("Neuron", "Astrocyte", "Oligodendrocyte"))
# print(subset_set)

# 8. Converting to bulk data

# Collapse to pseudobulk data (creates regular PhyloExpressionSet)
# bulk_set <- collapse(sc_phyex_set)
# print(bulk_set)

# 9. Transformations

# Apply log transformation to the single-cell data
# log_sc_set <- transform_counts(sc_phyex_set, log1p, "log1p")
# print(log_sc_set)

# Key features of ScPhyloExpressionSet:
# - Inherits from PhyloExpressionSet, so most bulk functions work
# - Automatically pseudobulks data by cell type for bulk-style analyses
# - Provides cell-level TAI calculations
# - Integrates with Seurat dimensional reductions
# - Specialized plotting functions for single-cell visualization
# - Statistical testing capabilities for comparing TAI across cell types
