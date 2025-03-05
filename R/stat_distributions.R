#' @import S7
Distribution <- new_class("Distribution",
      properties = list(
          name = new_required_property(
              class=class_character,
              name="name"
          ),
          pdf = new_required_property(
              class=class_function,
              name="pdf"
          ),
          cdf = new_required_property(
              class=class_function,
              name="cdf"
          ),
          quantile_function = new_required_property(
              class=class_function,
              name="quantile_function"
          ),
          fitting_function = new_required_property(
              class=class_function,
              name="fitting_function"
          ),
          param_names = new_required_property(
              class=class_character,
              name="param_names"
          )
          
      )
)



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


distributions <- list(normal = Distribution(name="normal",
                                            pdf=stats::dnorm,
                                            cdf=stats::pnorm,
                                            quantile_function=stats::qnorm,
                                            fitting_function=.fit_normal,
                                            param_names=c("mean", "sd")),
                      gamma  = Distribution(name="gamma",
                                            pdf=stats::dgamma,
                                            cdf=stats::pgamma,
                                            quantile_function=stats::qgamma,
                                            fitting_function=.fit_gamma,
                                            param_names=c("shape", "rate"))
)