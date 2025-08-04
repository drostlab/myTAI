#' @title Test Result S7 Class
#' @description S7 class for storing and manipulating statistical test results
#' from phylotranscriptomic conservation tests.
#' 
#' @param method_name Character string identifying the test method
#' @param test_stat Numeric test statistic value
#' @param fitting_dist Distribution object used for null hypothesis testing
#' @param params List of fitted distribution parameters
#' @param alternative Character string specifying alternative hypothesis ("two-sided", "less", "greater")
#' @param null_sample Numeric vector of null distribution samples
#' @param data_name Character string naming the dataset (optional)
#' @param p_label Character string for p-value label (default: "p_val")
#' 
#' @details
#' The TestResult class provides computed properties including:
#' - `p_value`: Computed p-value based on test statistic and fitted distribution
#' 
#' @examples
#' # Create a test result (typically done internally by conservation tests)
#' # result <- TestResult(method_name = "Test", test_stat = 1.5, ...)
#' 
#' @import S7
#' @export
TestResult <- new_class("TestResult",
    properties = list(
        ## CONSTRUCTOR PARAMETERS
        method_name = new_required_property(
            class=class_character,
            name="method_name"
        ),
        test_stat = new_required_property(
            class=class_numeric,
            name="test_stat"
        ),
        fitting_dist = new_required_property(
            class=Distribution,
            name="fitting_dist"
        ),
        params = new_required_property(
            class=class_list,
            name="params"
        ),
        alternative = new_options_property(
            class=class_character,
            options=c("two-sided", "less", "greater")
        ),
        null_sample = new_required_property(
            class = class_numeric,
            name="null_sample"
        ),
        data_name = new_property(
            class=class_character
        ),
        p_label = new_property(
            class=class_character,
            default="p_val"
        ),
        ## FIELDS & PROPERTIES
        p_value = new_property(
            class=class_double,
            getter= \(self) .get_p_value(cdf=self@fitting_dist@cdf,
                                        params=self@params,
                                        test_stat=self@test_stat,
                                        alternative=self@alternative)
        )

        ),
    validator = function(self) {
        # make sure params matches the fitting dist
        if (!setequal(names(self@params), self@fitting_dist@param_names))
            "the parameters provided do not match the fitting distribution used"
    }
    )

#' @title Calculate Confidence Intervals for Test Result
#' @description Calculate confidence intervals for the null distribution of a test result.
#' 
#' @param test_result A TestResult object
#' @param probs Numeric vector of probabilities for quantiles (default: c(0.025, 0.975))
#' 
#' @return Numeric vector of quantiles from the null distribution
#' 
#' @examples
#' # Calculate 95% confidence intervals
#' # ci <- conf_int(test_result, probs = c(0.025, 0.975))
#' 
#' @importFrom stats quantile
#' @keywords internal
conf_int <- function(test_result, 
                     probs=c(.025, .975)) {
    return(stats::quantile(test_result@null_sample, probs=probs))
}

#' @title Goodness of Fit Test
#' @description Perform a Kolmogorov-Smirnov test to assess goodness of fit
#' between the null sample and fitted distribution.
#' 
#' @param test_result A TestResult object
#' 
#' @return A ks.test result object
#' 
#' @details
#' This function tests whether the null sample follows the fitted distribution
#' using the Kolmogorov-Smirnov test. A significant result indicates poor fit.
#' 
#' @examples
#' # Test goodness of fit
#' # gof_result <- goodness_of_fit(test_result)
#' 
#' @keywords internal
goodness_of_fit <- function(test_result) {
    res <- do.call(ks.test,
                   c(list(x=test_result@null_sample,
                          y=test_result@fitting_dist@cdf),
                     test_result@params))
    return(res)
}


#' @import ggplot2
S7::method(plot, TestResult) <- function(test_result) {
    p <- ggplot(data.frame(x = test_result@null_sample), 
                aes(x = x)) +
        geom_histogram(
            aes(y = after_stat(density), fill = "Null Sample"), 
            bins = 100, 
            alpha = 0.7, 
            colour = "gray66") +
        geom_vline(
            aes(xintercept = test_result@test_stat, colour = "Test Statistic"), 
            linewidth = 1) +
        stat_function(fun = test_result@fitting_dist@pdf, args = test_result@params, aes(colour = "Fitted Null")) +
        scale_fill_manual(name = NULL, values = c("Null Sample" = "gray67")) +
        scale_colour_manual(name = NULL, values = c("Test Statistic" = "red", "Fitted Null" = "gray40")) +
        labs(x = "Score", y = "Density") +
        annotate("text",
                   x = test_result@test_stat + 0.05 * diff(range(test_result@null_sample)),
                   y = max(density(test_result@null_sample)$y) * 0.9,
                   label = exp_p(test_result@p_value),
                   parse=TRUE,
                   hjust = 1,
                   size = 3.5) +
        theme_minimal()
    
    return(p)
}

#' @title Print Method for TestResult
#' @description Print method for TestResult objects showing test summary information.
#' 
#' @param x A TestResult object
#' @param ... Additional arguments (ignored)
#' 
#' @return Invisibly returns the TestResult object
#' 
#' @examples
#' # Print a test result
#' # print(test_result)
#' 
#' @name print.TestResult
#' @keywords internal
S7::method(print, TestResult) <- function(x, ...) {
    cat("\n")
    cat("Statistical Test Result\n")
    cat("=======================\n")
    cat("Method:", x@method_name, "\n")
    cat("Test statistic:", x@test_stat, "\n")
    cat("P-value:", x@p_value, "\n")
    cat("Alternative hypothesis:", x@alternative, "\n")
    if (!is.null(x@data_name)) {
        cat("Data:", x@data_name, "\n")
    }
    cat("\n")
    invisible(x)
}


#' @title Plot Cullen-Frey Diagram for Distribution Assessment
#' @description Create a Cullen-Frey diagram to assess which distribution family
#' best fits the null sample data.
#' 
#' @param test_result A TestResult object
#' 
#' @return A Cullen-Frey plot from the fitdistrplus package
#' 
#' @details
#' The Cullen-Frey diagram plots skewness vs. kurtosis to help identify
#' appropriate distribution families for the null sample data.
#' 
#' @examples
#' # Create Cullen-Frey diagram
#' # cf_plot <- plot_cullen_frey(test_result)
#' 
#' @export
plot_cullen_frey <- function(test_result) {
    return(fitdistrplus::descdist(test_result@null_sample))
}



#' @title Conservation Test Result S7 Class
#' @description S7 class extending TestResult for conservation-specific test results,
#' including TXI profiles and null distributions.
#' 
#' @param method_name Character string specifying the statistical test method
#' @param test_stat Numeric value of the test statistic
#' @param fitting_dist Character string specifying the fitted distribution
#' @param params Named list of distribution parameters
#' @param alternative Character string specifying the alternative hypothesis
#' @param null_sample Numeric vector of null distribution values
#' @param data_name Character string describing the data
#' @param p_label Character string for p-value label
#' @param test_txi Numeric vector of observed TXI values
#' @param null_txis Matrix of null TXI distributions from permutations
#' @param modules Optional list of developmental modules used in the test
#' 
#' 
#' @details
#' ConservationTestResult extends TestResult with phylotranscriptomic-specific
#' information including the observed TXI profile and null TXI distributions
#' generated by permutation testing.
#' 
#' @examples
#' # Create a conservation test result (typically done internally)
#' # result <- ConservationTestResult(method_name = "Test", test_txi = c(1,2,3), ...)
#' 
#' @import S7
#' @export
ConservationTestResult <- new_class("ConservationTestResult",
    parent = TestResult,
    properties = list(
        test_txi = new_required_property(
            class = class_double,
            name="test_txi"
        ),
        null_txis = new_required_property(
            #class = class_matrix, # S7 doesn't support class_matrix yet
            name="null_txis"
        ),
        modules = new_property(
            class = class_list
        )
    )
    )

#' @title Plot Null TXI Sample Distribution
#' @description Create a plot showing the null TXI distribution sample compared 
#' to the observed test TXI values across developmental stages.
#' 
#' @param test_result A ConservationTestResult object containing null TXI distributions
#' 
#' @return A ggplot2 object showing null samples as gray lines and test TXI as colored line
#' 
#' @details
#' This function creates a visualization of the null hypothesis testing by plotting:
#' - Gray lines representing individual null TXI samples from permutations
#' - A horizontal line showing the mean of null distributions
#' - A colored line showing the observed test TXI values
#' 
#' The plot helps visualize how the observed TXI pattern compares to what would
#' be expected under the null hypothesis of no conservation signal.
#' 
#' @examples
#' # Plot null TXI sample distribution
#' # null_plot <- plot_null_txi_sample(conservation_test_result)
#' 
#' @import ggplot2 tibble tidyr dplyr RColorBrewer
#' @importFrom stats setNames
#' @importFrom tibble rowid_to_column
#' @export
plot_null_txi_sample <- function(test_result) {
    null_txis <- test_result@null_txis
    test_txis <- list("Test TXI" = test_result@test_txi)
    S <- ncol(null_txis)
    if (unique(sapply(test_txis, length)) != S)
        stop("The number of stages between the null_txis and the test_txis is inconsistent")
    
    null_df <- as_tibble(null_txis) |>
        rowid_to_column("Id") |>
        tidyr::pivot_longer(cols=-Id, names_to = "Stage", values_to = "TXI") |>
        add_column(Group = "Null Hypothesis Sample")
    
    test_df <- tibble(
        Id = rep(seq_along(test_txis), each = S),
        Group = rep(names(test_txis), each = S),
        Stage = rep(colnames(null_txis), times=length(test_txis)),
        TXI = unlist(test_txis)
    )
    
    avg <- mean(null_txis)
    
    colour_values <- stats::setNames(RColorBrewer::brewer.pal(length(unique(test_df$Group)), "Set2"), unique(test_df$Group))
    colour_values[1] <- "red"
    colour_values["Null Hypothesis Sample"] <- "gray67"
    
    p <- ggplot() +
        geom_line(data=null_df,
                  aes(x=factor(Stage, levels=unique(Stage)), 
                      y=TXI,  
                      group=interaction(Group, Id),
                      colour=factor(Group, levels=unique(Group))),
                  linewidth=0.1, alpha=0.1) + 
        geom_hline(yintercept=avg, colour="gray58", linewidth=0.7) +
        geom_line(data=test_df,
                  aes(x=factor(Stage, levels=unique(Stage)), 
                      y=TXI, 
                      colour=factor(Group, levels=unique(Group)), 
                      group=interaction(Group, Id)),
                  linewidth=1.6) +
        geom_line(data = null_df |> filter(Id == 1),
                  aes(x = factor(Stage, levels = unique(Stage)),
                      y = TXI,
                      colour = factor(Group, levels = unique(Group))),
                  linewidth = 0.5,
                  alpha = 0.7,
                  show.legend = TRUE) +
        scale_colour_manual(values=colour_values) +
        labs(x="Ontogeny",
             y="TAI",
             colour="Sample Type") + 
        scale_x_discrete(expand=c(0, 0)) +
        theme_minimal()
    
    return(p)
    
}


#' @title Calculate P-Value from Distribution
#' @description Internal function to calculate p-values from cumulative distribution functions
#' based on test statistics and alternative hypothesis specifications.
#' 
#' @param cdf Cumulative distribution function
#' @param test_stat Numeric test statistic value
#' @param params List of distribution parameters
#' @param alternative Character string specifying alternative hypothesis 
#'   ("two-sided", "less", "greater")
#' 
#' @return Numeric p-value
#' 
#' @details
#' This function calculates p-values using the appropriate tail(s) of the distribution:
#' - "greater": Uses upper tail (1 - CDF)
#' - "less": Uses lower tail (CDF)  
#' - "two-sided": Uses 2 * minimum of both tails
#' 
#' @examples
#' # Calculate p-value (internal use)
#' # pval <- .get_p_value(pnorm, 1.96, list(mean=0, sd=1), "two-sided")
#' 
#' @keywords internal
.get_p_value <- function(cdf, 
                         test_stat, 
                         params, 
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

