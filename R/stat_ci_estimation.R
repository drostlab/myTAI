
generate_bootstrapped_txis <- function(pTXI,
                                       count_matrix,
                                       sample_size) {
    
    message("generating bootstrapped txis")
    N <- nrow(pTXI)
    sampled_indices <- matrix(sample(N, size=N*sample_size, replace=TRUE), 
                              nrow=N, ncol=sample_size)
    selection_matrix <- Matrix::sparseMatrix(i = as.vector(sampled_indices),
                                             j = rep(1:sample_size, each = N),
                                             x = 1,
                                             dims= c(N, sample_size))
    
    # this is necessary to renormalise the counts so as to get the average right
    # it's a bit contrived, might be better to just rewrite this in rcpp...
    total_counts_per_bootstrap <- t(count_matrix) %*% selection_matrix |>
        as.matrix() |>
        t()
    total_counts <- colSums(count_matrix)
    
    bootstrap_matrix <- t(pTXI) %*% selection_matrix |> 
        as.matrix() |>
        t() |> # columns: conditions. rows: bootstraps 
        sweep(2, total_counts, "*")
    
    bootstrap_matrix <- bootstrap_matrix / total_counts_per_bootstrap
    
    colnames(bootstrap_matrix) <- colnames(pTXI)
    
    return(bootstrap_matrix)
}

memo_generate_bootstrapped_txis <- memoise::memoise(generate_bootstrapped_txis)





