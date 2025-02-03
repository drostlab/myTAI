

early_conservation.test <- function(phyex_set,
                                    modules,
                                    null_sample_size=10000,
                                    null_txis=NULL,
                                    plot_histogram=FALSE) {
    if(!setequal(c("early", "mid", "late"), names(modules)))
        stop('`modules` must have the structure: `list(early = ..., mid = ..., late = ...)`')
    
    # Null distribution: `ec_score` of null hypothesis TXI vectors,
    # obtained by permuting the age vector of the given `phyex_set`
    # These scores are fit to a normal distribution

    # 1. Simulate the TXI distribution under null hypothesis, via permutation
    # The user may supply existing parameters, which avoids expensive computation of permutations
    if (is.null(null_txis))
        null_txis <- generate_null_txi_sample(phyex_set, null_sample_size)
    
    # 2 - Compute the null distribution, using the `ec_score` of each TXI vector
    # as the test statistic
    null_sample <- apply(null_txis, 1, ec_score, modules=modules)
    
    # 3 - Fit a Normal distribution to the null sample of TXI variances
    null_dist_params <- .fit_normal(null_sample)
    
    # Compute test statistic: the `ec_score` of the TXI of the `phyex_set`
    txi <- TXI(phyex_set)
    ec_txi <- ec_score(txi, modules)
    
    # Get p value from observed variance and fitted null distribution
    pval <- stats::pnorm(ec_txi,
                         mean = null_dist_params$mean,
                         sd = null_dist_params$sd,
                         lower.tail = FALSE)
    
    if(plot_histogram)
        print(plot_ec_null_dist(null_sample, null_dist_params) +
                  ggplot2::geom_vline(xintercept=ec_txi, colour="red"))
    
    # Format the output to conform with the R stats hypothesis testing object `htest`
    
    mu_norm <- null_dist_params$mean
    names(mu_norm) <- "mean ec score"
    names(ec_txi) <- "ec score"
    rval <- list(statistic = ec_txi, parameters = null_dist_params, p.value = pval,
                 method = "Early Conservation Test", alternative="greater", null.value = mu_norm,
                 data.name=deparse(substitute(phyex_set)))
    class(rval) <- "htest"
    
    return(rval)
}

plot_ec_null_dist <- function(null_sample,
                              null_dist_params) {
    if (!all(c("mean", "sd") %in% names(null_dist_params)))
        stop("`null_dist_params` should be a list containing a `mean` and a `sd` parameter specifying a Normal distribution")
    
    # Plot gamma distribution
    p <- ggplot2::ggplot(data.frame(x = null_sample), ggplot2::aes(x=x)) +
        ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins=50, alpha=0.7) +
        ggplot2::stat_function(fun = dnorm, args = null_dist_params) +
        ggplot2::labs(x = "ec score", y = "Density") +
        ggplot2::theme_minimal()
    
    return(p)
}

ec_score <- function(txi, 
                     modules) {
    if(!setequal(c("early", "mid", "late"), names(modules)))
        stop('`modules` must have the structure: `list(early = ..., mid = ..., late = ...)`')
    D1 <- mean(txi[modules$mid]) - mean(txi[modules$early])
    D2 <- mean(txi[modules$late]) - mean(txi[modules$early])
    
    D_min <- min(D1, D2)
    
    return(D_min)
}

.fit_normal <- function(x) {
    params <- fitdistrplus::fitdist(x, "norm", method = "mme")
    return(list(mean=params$estimate[1], sd=params$estimate[2]))
}

