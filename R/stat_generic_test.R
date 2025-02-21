


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

#TODO return null samples and pdf used
generic_conservation_test <- function(phyex_set,
                                      test_name,
                                      scoring_function,
                                      fitting_dist,
                                      alternative) {
    
    
    # 1. Simulate the TXI distribution under null hypothesis, via permutation
    null_txis <- phyex_set@null_conservation_TXIs
    
    # 2. Compute the null distribution, using the appropriate scoring function
    null_sample <- apply(null_txis, 1, scoring_function)
    
    # 3. Compute test statistic using the same scoring function
    test_stat <- scoring_function(phyex_set@TXI)
    
    # 4. Estimate null distribution parameters
    if (fitting_dist == "gamma") {
        params <- .fit_gamma(null_sample)
        mu <- params$shape / params$rate
        cdf <- stats::pgamma
        pdf <- stats::dgamma
    }
    else if (fitting_dist == "normal") {
        params <- .fit_normal(null_sample)
        mu <- params$mean
        cdf <- stats::pnorm
        pdf <- stats::dnorm
    }
    else 
        stop("Currently supported fitting distributions are 'gamma' and 'normal'")
    
    # 5. Compute p value of the test
    pval <- get_p_value(cdf=cdf,
                        test_stat=test_stat,
                        params=params,
                        alternative=alternative)
    
    # Format the output to conform with the R stats hypothesis testing object `htest`
    # See https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/t.test.R
    # for an example
    names(mu) <- "mean score"
    names(test_stat) <- "score"
    rval <- list(statistic = test_stat, parameters = params, p.value = pval,
                 method = test_name, alternative=alternative, null.value = mu,
                 data.name=deparse(substitute(phylo_set)))
    
    #TODO: compute confidence interval for test
    
    class(rval) <- "htest"
    
    return(rval)
}