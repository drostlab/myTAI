
#' @title Destroy Phylotranscriptomic Pattern Using GATAI
#' @description Apply the GATAI algorithm
#' to identify and remove genes that contribute to phylotranscriptomic patterns.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param num_runs Number of GATAI runs to perform (default: 20)
#' @param runs_threshold Threshold for gene removal consistency across runs (default: 0.5)
#' @param analysis_dir Directory to store GATAI analysis results (default: NULL)
#' @param max_generations: Integer. Maximum number of generations (iterations) for the genetic algorithm (default 10000).
#' @param seed Random seed for reproducibility (default: 1234)
#' @param ... Additional arguments passed to gataiR::gatai
#' 
#' @return A list containing GATAI results including identified genes that contribute to the pattern
#' 
#' @details
#' This function requires the gataiR package to be installed. GATAI systematically removes genes
#' that contribute to phylotranscriptomic patterns by iteratively testing gene removal and
#' evaluating the impact on the overall transcriptomic signature.
#' 
#' @examples
#' # Apply GATAI to identify pattern-contributing genes
#' # gatai_result <- destroy_pattern(phyex_set, num_runs = 10, runs_threshold = 0.6)
#' 
#' @author Filipa Martins Costa
#' @export
destroy_pattern <- function(phyex_set, 
                            num_runs = 20,
                            runs_threshold = 0.5,
                            analysis_dir = NULL,
                            max_generations = 10000,
                            seed = 1234,
                            ...) {
    if (!requireNamespace("gataiR", quietly = TRUE)) {
        stop("Package 'gataiR' must be installed to use this function.")
    }
    
    res <- gataiR::gatai(phyex_set@data_collapsed, 
                         num_runs = num_runs,
                         runs_threshold = runs_threshold, 
                         max_generations = max_generations,
                         seed = seed,
                         ...)
    res <- list(removed_genes=res$common_removed_genes, runs=res$genes_list)
    

    
    return(res)
}


#' @title Plot Comprehensive GATAI Results
#' @description Create a suite of plots summarizing the effects of GATAI gene removal on phylotranscriptomic patterns.
#'
#' @param phyex_set A PhyloExpressionSet object containing the original gene expression data.
#' @param gatai_result Result list from \code{destroy_pattern()}, containing GATAI analysis output.
#' @param conservation_test Function for conservation test (default: \code{flatline_test}).
#' @param runs_threshold Threshold for gene removal consistency across runs (default: 0.5).
#' @param signature_plot_type Type of signature plot: "separate" for individual plots, "combined" for overlay (default: both options).
#' @param n_random_genes Number of random genes to remove for comparison (default: same as GATAI removed).
#'
#' @return A named list of ggplot/cowplot objects and results:
#'   \item{signature_plots}{Signature plots before/after GATAI and random removal}
#'   \item{heatmap_plot}{Heatmap of GATAI-removed genes}
#'   \item{profiles_plot}{Gene expression profiles of GATAI-removed genes}
#'   \item{profiles_plot_facet}{Faceted gene profiles by strata}
#'   \item{gene_space_plot}{Gene space plot of GATAI-removed genes}
#'   \item{mean_var_plot}{Mean-variance plot highlighting GATAI-removed genes}
#'   \item{strata_plot}{Phylostrata distribution plot (log obs/exp) for GATAI-removed genes}
#'   \item{null_dist_plot}{Null distribution plot with test statistics and p-values}
#'   \item{convergence_plots}{GATAI convergence plots (if available)}
#'
#' @details
#' This function provides a comprehensive visualization of the impact of GATAI gene removal,
#' including transcriptomic signature plots, gene expression profiles, heatmaps, mean-variance relationships,
#' phylostrata distributions, conservation test comparisons, and convergence diagnostics.
#'
#' @examples
#' # Run GATAI analysis
#' # gatai_result <- destroy_pattern(phyex_set, num_runs = 20)
#' # Plot results
#' # plots <- plot_gatai_results(phyex_set, gatai_result)
#' # Print signature plots
#' # print(plots$signature_plots)
#'
#' @author Filipa Martins Costa, Stefan Manolache
#' @export
plot_gatai_results <- function(phyex_set,
                               gatai_result,
                               conservation_test = flatline_test,
                               runs_threshold = 0.5,
                               signature_plot_type = c("separate", "combined"),
                               n_random_genes = length(gatai_result$removed_genes)) {
    plot_type <- match.arg(signature_plot_type)
    removed_genes <- gatai_result$removed_genes   

    # 1. Plot TAI signature before/after GATAI
    gatai_set <- remove_genes(phyex_set, removed_genes, new_name = paste(phyex_set@name, "- GATAI removed"))
    random_removed_genes <- sample(phyex_set@gene_ids, n_random_genes)
    random_set <- remove_genes(phyex_set, random_removed_genes, new_name = paste(phyex_set@name, "- Random removed"))

    if (plot_type == "separate") {
        P1 <- plot_signature(phyex_set, show_p_val = TRUE, conservation_test = conservation_test, colour="blue") +
            ggtitle(paste("Original:", phyex_set@num_genes, "genes"))
        
        P2 <- plot_signature(gatai_set, show_p_val = TRUE, conservation_test = conservation_test, colour="red") +
            ggtitle(paste("GATAI removed:", length(removed_genes), "genes"))
        
        P3 <- plot_signature(random_set, show_p_val = TRUE, conservation_test = conservation_test, colour="darkgray") +
            ggtitle(paste("Random removed:", n_random_genes, "genes"))
        # Align y-axis scales and combine plots
        y_limits <- range(c(phyex_set@TXI, gatai_set@TXI, random_set@TXI)) + c(-0.05, 0.05)
        plots <- list(P1, P2, P3)
        plots <- lapply(plots, \(p) p + ylim(y_limits))
        signature_plots <- cowplot::plot_grid(plotlist = plots, ncol = 3)
    }
    else {
        signature_plots <- plot_signature_multiple(c(phyex_set, gatai_set, random_set), 
                                                   show_p_val = TRUE, conservation_test = conservation_test, colours=c("blue", "red", "darkgray")) +
            ggtitle(paste("Signature Comparison: Original vs GATAI-removed (", 
                          length(removed_genes), "genes) vs Random-removed (", n_random_genes, "genes)"))
    }


    # 2. Plot removed gene profiles and heatmap
    heatmap_plot <- plot_gene_heatmap(phyex_set, genes=removed_genes, cluster_rows=TRUE, show_gene_ids=TRUE) +
        ggtitle(paste("Expression Heatmap: GATAI-removed genes (", length(removed_genes), "genes)"))
    profiles_plot <- plot_gene_profiles(phyex_set, genes=removed_genes, max_genes=200, transformation="none") +
        ggtitle(paste("Gene Expression Profiles: GATAI-removed genes (max 200 shown)"))
    profiles_plot_facet <- plot_gene_profiles(phyex_set, genes=removed_genes, max_genes=1000, facet_by_strata = TRUE, transformation="none") +
        ggtitle(paste("Gene Expression Profiles by Phylostrata: GATAI-removed genes")) +
        theme(strip.text = element_text(size = 6))

    gene_space_plot <- plot_gene_space(phyex_set, genes=removed_genes, colour_by="strata") +
        ggtitle(paste("Gene Space Plot: GATAI-removed genes (", length(removed_genes), "genes)"))

    # 3. Mean-variance plot highlighting GATAI-removed genes
    mean_var_plot <- plot_mean_var(phyex_set, highlight_genes = removed_genes) +
        ggtitle(paste("Mean-Variance Plot with GATAI-removed genes highlighted (", length(removed_genes), "genes)"))
    
    # 4. Distribution of gatai genes by age with log obs/exp ratios
    strata_plot <- plot_distribution_strata(phyex_set@strata,
                                            selected_gene_ids = removed_genes,
                                            as_log_obs_exp = TRUE) +
        ggtitle("Phylostrata Distribution: GATAI-removed genes (Log Obs/Exp)")
    
    # 5. Compare p values
    original_test <- conservation_test(phyex_set, plot_result = FALSE)
    gatai_test <- conservation_test(gatai_set, plot_result = FALSE)
    
    null_sample <- original_test@null_sample
    test_stats <- data.frame(
        label = c("Original", "GATAI removed"),
        stat = c(original_test@test_stat, gatai_test@test_stat),
        p_value = c(original_test@p_value, gatai_test@p_value),
        color = c("blue", "red")
    )

    null_dist_plot <- ggplot(data.frame(x = null_sample), aes(x = x)) +
        geom_histogram(aes(y = after_stat(density)), bins = 100, fill = "gray67", alpha = 0.7, colour = "gray66") +
        stat_function(fun = original_test@fitting_dist@pdf, args = original_test@params, colour = "gray40") +
        geom_vline(data = test_stats, aes(xintercept = stat, colour = label), linewidth = 1) +
        scale_colour_manual(name = NULL, values = c("Original" = "blue", "GATAI removed" = "red")) +
        labs(x = "Score", y = "Density",
             title = "Conservation Test Comparison",
             subtitle = paste("Test:", original_test@method_name)) +
        annotate("text",
                 x = test_stats$stat - 0.05 * diff(range(null_sample)),
                 y = max(density(null_sample)$y) * c(0.9, 0.8),
                 label = sapply(test_stats$p_value, function(p) formatC(p, format = "e", digits = 2)),
                 colour = test_stats$color,
                 hjust = 1,
                 size = 3.5) +
        theme_minimal()

    

    # 6. Convergence plots
    convergence_plots <- full_gatai_convergence_plot(phyex_set, gatai_result$runs, p=runs_threshold) +
        ggtitle(paste("GATAI Convergence Analysis (threshold:", runs_threshold, ")")) +
        theme(aspect.ratio = 12/8,
              text = element_text(size = 7),
              axis.text = element_text(size = 6),
              plot.title = element_text(size = 9),
              legend.text = element_text(size = 6),
              legend.title = element_text(size = 7))

    result_list <- list(
        signature_plots = signature_plots,
        heatmap_plot = heatmap_plot,
        profiles_plot = profiles_plot,
        profiles_plot_facet = profiles_plot_facet,
        gene_space_plot = gene_space_plot,
        mean_var_plot = mean_var_plot,
        strata_plot = strata_plot,
        null_dist_plot = null_dist_plot,
        convergence_plots = convergence_plots
    )
    return(result_list)
}

#' @title Save GATAI Analysis Results to PDF
#' @description Save removed gene IDs and all GATAI analysis plots to a PDF file.
#'
#' @param analysis_dir Directory to save the PDF file.
#' @param phyex_set A PhyloExpressionSet object containing the original gene expression data.
#' @param gatai_result Result list from \code{destroy_pattern()}, containing GATAI analysis output.
#' @param prefix Optional prefix for the PDF filename (default: "GATAI_analysis").
#' @param ... Additional arguments passed to \code{plot_gatai_results()}.
#'
#' @return Invisibly returns the path to the saved PDF.
#'
#' @examples
#' # Save results after running destroy_pattern
#' # save_gatai_results_pdf("results/", phyex_set, gatai_result)
#'
#' @export
save_gatai_results_pdf <- function(analysis_dir = "gatai_analysis",
                                   phyex_set,
                                   gatai_result,
                                   prefix = "GATAI_analysis",
                                   ...) {
    if (!dir.exists(analysis_dir)) {
        dir.create(analysis_dir, recursive = TRUE)
    }
    pdf_path <- file.path(analysis_dir, paste0(prefix, ".pdf"))
    # Generate plots using plot_gatai_results
    plots <- plot_gatai_results(phyex_set, gatai_result, ...)
    removed_genes <- gatai_result$removed_genes
    grDevices::pdf(pdf_path, width = 10, height = 8)
    
    # Page 1: Removed gene IDs - handle long gene lists better
    plot.new()
    title("GATAI Removed Gene IDs")
    text(0, 0.95, paste("Total removed genes:", length(removed_genes)), adj = c(0, 1), cex = 1.2)
    text(0, 0.92, "Gene IDs:", adj = c(0, 1), cex = 1)
    
    # Better handling of long gene lists
    gene_text <- paste(removed_genes, collapse = ", ")
    wrapped <- strwrap(gene_text, width = 100)  # Reduced width for better fit
    
    # Calculate how many lines can fit on page
    max_lines <- floor((0.90 - 0.1) / 0.02)  # Leave space at bottom and adjust for new spacing
    
    if (length(wrapped) > max_lines) {
        # Split across multiple pages
        pages_needed <- ceiling(length(wrapped) / max_lines)
        
        for (page in 1:pages_needed) {
            if (page > 1) {
                plot.new()
                title(paste("GATAI Removed Gene IDs (continued, page", page, "of", pages_needed, ")"))
            }
            
            start_line <- (page - 1) * max_lines + 1
            end_line <- min(page * max_lines, length(wrapped))
            
            for (i in start_line:end_line) {
                line_pos <- i - start_line + 1
                text(0, 0.90 - 0.02 * line_pos, wrapped[i], adj = c(0, 1), cex = 0.8)
            }
        }
    } else {
        # All fits on one page
        for (i in seq_along(wrapped)) {
            text(0, 0.90 - 0.02 * i, wrapped[i], adj = c(0, 1), cex = 0.8)
        }
    }
    
    # Remaining pages: Plots without titles
    for (nm in names(plots)) {
        try({
            print(plots[[nm]])
            # Remove the title() call to avoid showing plot names
        }, silent = TRUE)
    }
    
    grDevices::dev.off()
    message("Saved GATAI analysis PDF to: ", pdf_path)
    invisible(pdf_path)
}