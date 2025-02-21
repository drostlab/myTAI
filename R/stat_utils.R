

#' @import ggplot2
#' TODO plot confidence interval limits?
#' TODO plot from htest object
plot_null_dist <- function(null_sample,
                           pdf,
                           params) {
    
    p <- ggplot(data.frame(x = null_sample), 
                aes(x=x)) +
        geom_histogram(
            aes(y = ggplot2::after_stat(density)), 
            bins=50, 
            alpha=0.7, 
            fill="gray67", 
            colour="gray66") +
        stat_function(fun = pdf, args = params, colour="gray40") +
        labs(x = "Score", y = "Density") +
        theme_minimal()
    
    return(p)
}

get_p_value <- function(cdf, 
                        test_stat, 
                        params = list(), 
                        alternative = c("two-sided", "less", "greater")) {
    alternative <- match.arg(alternative)
    
    cdf_value_low <- do.call(cdf, c(list(test_stat, lower.tail=TRUE), params))
    cdf_value_high <- do.call(cdf, c(list(test_stat, lower.tail=FALSE), params))
    
    pval <- switch(
        alternative,
        "greater" = cdf_value_high,
        "less" = cdf_value_low,
        "two-sided" = 2 * min(cdf_value_low, cdf_value_high),
    )
    return(pval)
}

.fit_normal <- function(x) {
    params <- fitdistrplus::fitdist(x, "norm", method = "mme")
    return(list(mean=params$estimate[1], sd=params$estimate[2]))
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