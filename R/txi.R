

TXI <- function(phyex_set) {
    age_vector <- as.vector(phyex_set[[1]])
    count_matrix <- as.matrix(phyex_set[3:ncol(phyex_set)])

    txi <- cpp_TAI(count_matrix, age_vector)
    
    return(txi)
}

TAI <- TXI
TDI <- TXI
