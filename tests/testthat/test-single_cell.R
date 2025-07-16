test_that("get_sc_TAI works with mock data", {
    # Skip if Seurat is not available
    skip_if_not_installed("Seurat")
    
    # Create mock Seurat object
    mock_counts <- matrix(
        rpois(1000, lambda = 10),
        nrow = 100,
        ncol = 10,
        dimnames = list(
            paste0("Gene", 1:100),
            paste0("Cell", 1:10)
        )
    )
    
    # Create mock phylomap
    phylomap <- data.frame(
        Stratum = sample(1:5, 100, replace = TRUE),
        GeneID = paste0("Gene", 1:100)
    )
    
    # Create minimal Seurat object
    seurat_obj <- Seurat::CreateSeuratObject(counts = mock_counts)
    
    # Test get_sc_TAI
    sc_tai <- get_sc_TAI(seurat_obj, phylomap)
    
    expect_true(is.numeric(sc_tai))
    expect_equal(length(sc_tai), ncol(mock_counts))
    expect_true(all(is.finite(sc_tai)))
    expect_true(all(sc_tai >= 0))
})

test_that("get_sc_TAI handles different slots", {
    # Skip if Seurat is not available
    skip_if_not_installed("Seurat")
    
    # Create mock data
    mock_counts <- matrix(
        rpois(1000, lambda = 10),
        nrow = 100,
        ncol = 10,
        dimnames = list(
            paste0("Gene", 1:100),
            paste0("Cell", 1:10)
        )
    )
    
    phylomap <- data.frame(
        Stratum = sample(1:5, 100, replace = TRUE),
        GeneID = paste0("Gene", 1:100)
    )
    
    seurat_obj <- Seurat::CreateSeuratObject(counts = mock_counts)
    
    # Test with different slots
    sc_tai_data <- get_sc_TAI(seurat_obj, phylomap, slot = "data")
    sc_tai_counts <- get_sc_TAI(seurat_obj, phylomap, slot = "counts")
    
    expect_true(is.numeric(sc_tai_data))
    expect_true(is.numeric(sc_tai_counts))
    expect_equal(length(sc_tai_data), length(sc_tai_counts))
})

test_that("get_sc_TAI handles missing genes", {
    # Skip if Seurat is not available
    skip_if_not_installed("Seurat")
    
    # Create mock data
    mock_counts <- matrix(
        rpois(1000, lambda = 10),
        nrow = 100,
        ncol = 10,
        dimnames = list(
            paste0("Gene", 1:100),
            paste0("Cell", 1:10)
        )
    )
    
    # Create phylomap with some genes not in the counts
    phylomap <- data.frame(
        Stratum = sample(1:5, 150, replace = TRUE),
        GeneID = paste0("Gene", 1:150)  # 50 genes not in counts
    )
    
    seurat_obj <- Seurat::CreateSeuratObject(counts = mock_counts)
    
    # Should work despite missing genes
    sc_tai <- get_sc_TAI(seurat_obj, phylomap)
    expect_true(is.numeric(sc_tai))
    expect_equal(length(sc_tai), ncol(mock_counts))
})
