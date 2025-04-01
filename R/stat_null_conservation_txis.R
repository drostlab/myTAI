


generate_conservation_txis <- function(strata_vector,
                                       count_matrix,
                                       sample_size) {
    
    permuted_txi_matrix <- cpp_bootMatrix(count_matrix,
                                          strata_vector,
                                          sample_size)
    
    colnames(permuted_txi_matrix) <- colnames(count_matrix)
    return(permuted_txi_matrix)
}
memo_generate_conservation_txis <- memoise::memoise(generate_conservation_txis)