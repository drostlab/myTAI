


flt_p_val <- function(phyex_set,
                      null_dist) {
    is.ExpressionSet(phyex_set)
    
    age_vector <- as.vector(phyex_set[[1]])
    count_matrix <- as.matrix(phyex_set[3:ncol(phyex_set)])
    # compute TAI variance
    tai <- cpp_TAI(count_matrix, age_vector)
    var <- stats::var(tai)
    
    # get p value from observed variance and fitted null distribution
    p_val <- stats::pgamma(var,
                          shape = null_dist$shape,
                          rate = null_dist$rate,
                          lower.tail = FALSE)
    return(p_val)
}

flt_null_dist <- function(phyex_set, num_samples = 10000) {
    # a sample of TAI variances
    sample <- flt_null_dist_sample(phyex_set, num_samples)
    
    # fit a gamma distribution based on the variances
    params <- GetGamma(sample)
    
    return(list(shape=params[[1]], rate=params[[2]]))
}

flt_null_dist_sample <- function(phyex_set,
                                  num_samples) {
    is.ExpressionSet(phyex_set)
    
    age_vector <- as.vector(phyex_set[[1]])
    count_matrix <- as.matrix(phyex_set[3:ncol(phyex_set)])
    
    # columns: stages, rows: permuted instances, entries: TAI value
    permuted_tai_matrix <- cpp_bootMatrix(count_matrix,
                                          age_vector,
                                          num_samples)
    
    permuted_vars <- apply(permuted_tai_matrix, 1, stats::var)
    
    
    return(permuted_vars)
}