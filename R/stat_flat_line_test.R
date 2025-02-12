

diagnose_test_robustness <- function(phyex_set, test) {
    num_runs = 10
    sample_sizes = c(1000, 10000)
    
    f <- function(size) {
        return(test(phyex_set, null_sample_size=size))
    }
    # make parallel
    res_vec <- purrr::map(rep(sample_sizes, each=10), f) |>
        purrr::map("p.value")
        
    return(res_vec)
}

# (Stefan) - I used the *.test naming convention, following R stats
stat_flat_line_test <- function(phylo_set,
                           null_sample_size=10000,
                           null_txis=NULL,
                           plot_histogram=FALSE) {
    
    # Null distribution: variances of null hypothesis TXI vectors,
    # obtained by permuting the age vector of the given `phyex_set`
    # These variances are fit to a gamma distribution
    
    # 1. Simulate the TXI distribution under null hypothesis, via permutation
    # The user may supply an existing null txi sample
    if (is.null(null_txis))
        null_txis <- generate_null_txi_sample(phylo_set, null_sample_size)
    ##TODO check null_txis is valid
    
    # 2. Compute the null distribution, using the variance of each TXI vector
    null_sample <- apply(null_txis, 1, stats::var)
    
    # 3. Fit a Gamma distribution to the null sample of TXI variances
    null_dist_params <- .fit_gamma(null_sample)
    
    # Compute test statistic: the variance of the TXI of the `phyex_set`
    txi <- TXI(phylo_set)
    var_txi <- stats::var(txi)
    
    # Plot distribution
    if(plot_histogram) 
        print(plot_flt_null_dist(null_sample, null_dist_params) +
              ggplot2::geom_vline(xintercept=var_txi, colour="red"))
    
    # Get p value from observed variance and fitted null distribution
    pval <- stats::pgamma(var_txi,
                          shape = null_dist_params$shape,
                          rate = null_dist_params$rate,
                          lower.tail = FALSE)
    
    # Format the output to conform with the R stats hypothesis testing object `htest`
    # See https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/t.test.R
    # for an example
    
    mu_gamma <- null_dist_params$shape / null_dist_params$rate
    names(mu_gamma) <- "mean var"
    names(var_txi) <- "var"
    rval <- list(statistic = var_txi, parameters = null_dist_params, p.value = pval,
                 method = "Flat Line Test", alternative="greater", null.value = mu_gamma,
                 data.name=deparse(substitute(phylo_set)))
    class(rval) <- "htest"
    
    return(rval)
}


plot_flt_null_dist <- function(null_sample,
                               null_dist_params) {
    if (!all(c("shape", "rate") %in% names(null_dist_params)))
        stop("`null_dist_params` should be a list containing a `shape` and a `rate` parameter specifying a Gamma distribution")
    
    # Plot gamma distribution
    p <- ggplot2::ggplot(data.frame(x = null_sample), ggplot2::aes(x=x)) +
        ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins=50, alpha=0.7, fill="gray67", colour="gray66") +
        ggplot2::stat_function(fun = dgamma, args = null_dist_params, colour="gray40") +
        ggplot2::labs(x = "Variance", y = "Density") +
        ggplot2::theme_minimal()
    
    return(p)
}

.fit_gamma <- function(x) {
    iterations = 200
    max_cut = 0.25
    step = max_cut/iterations 
    sorted_vars = sort(x, decreasing = TRUE)
    max_p_fit_v = 0
    max_p_i = 0
    for (i in 2:iterations) {
        # to avoid indexing from zero
        # Filtered variances
        filtered_vars <-
            sorted_vars[round(length(x) * i * step):length(x)]
        
        # Estimate parameters using method of moments
        gamma_fit <-
            fitdistrplus::fitdist(filtered_vars, "gamma", method = "mme")
        shape <- gamma_fit$estimate[1]
        rate <- gamma_fit$estimate[2]
        # Perform Kolmogorov-Smirnov test
        suppressWarnings(ks_result <-
                             stats::ks.test(filtered_vars,
                                            "pgamma",
                                            shape = shape,
                                            rate = rate))
        if (ks_result$p.value > max_p_fit_v) {
            max_p_i = i
            max_p_fit_v = ks_result$p.value
        }
    }
    if (max_p_i == 0) {
        gamma_fit <- fitdistrplus::fitdist(x, "gamma", method = "mme")
        return(list(shape=gamma_fit$estimate[1], 
                    rate=gamma_fit$estimate[2]))
    }
    b_shape = 0
    b_rate = 0
    ks_best = NULL
    max_p_fit_v = 0
    for (i in -10:10) {
        # Filtered variances
        lb = round(length(x) * (max_p_i * step + i * step / 10))
        filtered_vars <-
            sorted_vars[lb:length(x)]
        
        # Estimate parameters using method of moments
        gamma_fit <-
            fitdistrplus::fitdist(filtered_vars, "gamma", method = "mme")
        shape <- gamma_fit$estimate[1]
        rate <- gamma_fit$estimate[2]
        # Perform Kolmogorov-Smirnov test
        suppressWarnings(ks_result <-
                             stats::ks.test(filtered_vars,
                                            "pgamma",
                                            shape = shape,
                                            rate = rate))
        if (ks_result$p.value > max_p_fit_v) {
            max_p_fit_v = ks_result$p.value
            b_shape = shape
            b_rate = rate
            ks_best = ks_result
        }
    }
    return(list(shape=b_shape, rate=b_rate))
}

