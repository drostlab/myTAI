generic_conservation_test <- function(phyex_set,
                                      test_name,
                                      scoring_function,
                                      fitting_dist,
                                      alternative = c("two-sided", "greater", "less"),
                                      custom_null_txis = NULL
                                      ) {
    # check arguments
    S7::check_is_S7(phyex_set, PhyloExpressionSet)
    stopifnot(is.function(scoring_function))
    alternative <- match.arg(alternative)
    
    # 1. Simulate the TXI distribution under null hypothesis, via permutation
    if (!is.null(custom_null_txis))
        null_txis <- custom_null_txis
    else
        null_txis <- phyex_set@null_conservation_txis
    
    # 2. Compute the null distribution, using the appropriate scoring function
    null_sample <- apply(null_txis, 1, scoring_function)
    
    # 3. Compute test statistic using the same scoring function
    test_stat <- scoring_function(phyex_set@TXI)
    
    # 4. Estimate null distribution parameters
    params <- fitting_dist@fitting_function(null_sample)
    
    res <- ConservationTestResult(method_name=test_name,
                                  test_stat=test_stat,
                                  fitting_dist=fitting_dist,
                                  params=params,
                                  alternative=alternative,
                                  null_sample=null_sample,
                                  data_name=deparse(substitute(phylo_set)),
                                  null_txis=null_txis,
                                  test_txi=phyex_set@TXI
                                  )
    
   return(res)
}

#' @import ggplot2
diagnose_test_robustness <- function(test, 
                                     phyex_set,
                                     sample_sizes=c(500, 1000, 5000, 10000),
                                     num_reps=5,
                                     ...) {

    f <- function(size) {
        null_txis <- .generate_conservation_txis(phyex_set, sample_size=size)
        return(test(phyex_set, custom_null_txis=null_txis, ...))
    }
    
    res_vec <- purrr::map(rep(sample_sizes, each=num_reps), f) |>
        purrr::map_dbl(~ .x@p_value)
    
    df <- data.frame(pval=res_vec, sample_size=rep(sample_sizes, each=num_reps))
    
    p <- ggplot(df, aes(x=factor(sample_size), y=pval, color="Independent run")) +
        geom_jitter(height=0, width=0.01) +
        scale_x_discrete(breaks=sample_sizes) +
        scale_y_continuous(transform='log10') +
        labs(x="Null sample size", 
             y="p value",
             color="Legend",
             title="Robustness of test p-values against different null distribution sample sizes") +
        scale_color_manual(values="black") +
        theme_minimal()
    
    return(p)
}
