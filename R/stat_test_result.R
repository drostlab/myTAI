#' @import S7
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
        ## FIELDS & PROPERTIES
        p_value = new_cached_property(
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

conf_int <- function(test_result, 
                     probs=c(.025, .975)) {
    return(quantile(test_result@null_sample, probs=probs))
}

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
                aes(x=x)) +
        geom_histogram(
            aes(y = after_stat(density)), 
            bins=100, 
            alpha=0.7, 
            fill="gray67", 
            colour="gray66") +
        geom_vline(
            xintercept=test_result@test_stat,
            colour="red") +
        stat_function(fun = test_result@fitting_dist@pdf, args = test_result@params, colour="gray40") +
        labs(x = "Score", y = "Density") +
        theme_minimal()
    
    return(p)
}

# S7::S4_register(TestResult)
# S7::method(print, TestResult) <- function(x) {
#     return("a")
# }


plot_cullen_frey <- function(test_result) {
    return(fitdistrplus::descdist(test_result@null_sample))
}



#' @import S7
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

# plot samples against true txi.
#' @import ggplot2 tibble
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
        Id = rep(1:length(test_txis), each = S),
        Group = rep(names(test_txis), each = S),
        Stage = rep(colnames(null_txis), times=length(test_txis)),
        TXI = unlist(test_txis)
    )
    
    avg <- mean(null_txis)
    
    colour_values <- setNames(RColorBrewer::brewer.pal(length(unique(test_df$Group)), "Set2"), unique(test_df$Group))
    colour_values["Null Hypothesis Sample"] <- "paleturquoise3"
    
    p <- ggplot(test_df,
                aes(x=factor(Stage, levels=unique(Stage)), 
                    y=TXI, 
                    colour=factor(Group, levels=unique(Group)), 
                    group=interaction(Group, Id))) +
        geom_line(linewidth=1.6) + 
        geom_line(data=null_df,
                  aes(x=factor(Stage, levels=unique(Stage)), 
                      y=TXI,  
                      group=interaction(Group, Id),
                      colour=factor(Group, levels=unique(Group))),
                  linewidth=0.1, alpha=0.1) + 
        geom_hline(yintercept=avg, colour="paleturquoise4", linewidth=0.7) +
        scale_colour_manual(values=colour_values) +
        labs(x="Ontogeny",
             y="TAI",
             colour="Sample Type") + 
        scale_x_discrete(expand=c(0, 0)) +
        theme_minimal()
    
    return(p)
    
}


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

