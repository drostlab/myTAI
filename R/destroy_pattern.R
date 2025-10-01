#' @title Destroy Phylotranscriptomic Pattern Using GATAI
#' @description Apply the GATAI algorithm
#' to identify and remove genes that contribute to phylotranscriptomic patterns.
#' 
#' @param phyex_set A PhyloExpressionSet object (bulk or single cell, the latter which will get pseudo-bulked)
#' @param num_runs Number of GATAI runs to perform (default: 20)
#' @param runs_threshold Threshold for gene removal consistency across runs (default: 0.5)
#' @param analysis_dir Directory to store GATAI analysis results (default: NULL)
#' @param plot_results Whether to plot the results. If analysis dir is given, this will be ignored.
#' @param max_generations Integer. Maximum number of generations (iterations) for the genetic algorithm (default 10000).
#' @param seed Random seed for reproducibility (default: 1234)
#' @param extended_analysis Whether to show the multiple runs and convergence plots (default: FALSE)
#' @param min_pval Minimum p-value for which the pattern is considered destroyed (default: 0.05).
#' @param always_return_genes Whether to return genes even when the pattern is not destroyed (default: FALSE).
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
#' @import patchwork
#' @author Filipa Martins Costa
#' @export
destroy_pattern <- function(phyex_set, 
                            num_runs = 20,
                            runs_threshold = 0.5,
                            analysis_dir = NULL,
                            plot_results = TRUE,
                            max_generations = 10000,
                            seed = 1234,
                            extended_analysis = FALSE,
                            min_pval = 0.05,
                            always_return_genes = FALSE,
                            ...) {
    if (!requireNamespace("gataiR", quietly = TRUE)) {
        stop("Package 'gataiR' must be installed to use this function.")
    }
    check_PhyloExpressionSet(phyex_set)
    
    gatai_res <- gataiR::gatai(as_data_frame(collapse(phyex_set)), 
                         num_runs = num_runs,
                         runs_threshold = runs_threshold, 
                         max_generations = max_generations,
                         seed = seed,
                         ...)
    res <- list(removed_genes=gatai_res$common_removed_genes)
    if (extended_analysis)
        res$runs <- gatai_res$genes_list

    if (length(res$removed_genes) == 0) {
        message("GATAI has failed to detect any genes. Try to increase the number of generations.")
        return(res)
    }

    
    # compute p value of destroyed set
    pval <- phyex_set |> remove_genes(res$removed_genes) |> stat_flatline_test(plot_result=FALSE) |> (\(t) t@p_value)()
    if (pval < min_pval) {
        message("GATAI has failed to destroy the pattern. Try specifying different parameters.")
        if (!always_return_genes) {
            res$removed_genes <- character(0)
            return(res)
        }
        else {
            message("Returning the genes anyway...")
        }
    }
    

    
    # If analysis_dir is provided, save results to PDF and genes.txt
    if (!is.null(analysis_dir)) {
        message(sprintf("Saving results to %s", analysis_dir))
        # Save genes.txt
        if (!dir.exists(analysis_dir)) {
            dir.create(analysis_dir, recursive = TRUE)
        }
        genes_path <- file.path(analysis_dir, "genes.txt")
        writeLines(res$removed_genes, genes_path)
        
        # If extended_analysis is TRUE and runs are available, save individual run gene files
        if (extended_analysis && "runs" %in% names(res)) {
            runs_dir <- file.path(analysis_dir, "individual_runs")
            if (!dir.exists(runs_dir)) {
                dir.create(runs_dir, recursive = TRUE)
            }
            
            for (i in seq_along(res$runs)) {
                run_genes_path <- file.path(runs_dir, paste0("run_", i, "_genes.txt"))
                writeLines(res$runs[[i]], run_genes_path)
            }
            message(sprintf("Saved individual run gene files to %s", runs_dir))
        }
        
        save_gatai_results_pdf(phyex_set = phyex_set,
                               gatai_result = res,
                               analysis_dir = analysis_dir,
                               runs_threshold = runs_threshold)
    }
    else {
        message("No analysis directory name was provided. To save the results to file, pass `analysis_dir`")
        if (plot_results)
            res$plots <- plot_gatai_results(phyex_set, res, runs_threshold=runs_threshold)
    }
    
    return(res)
}


#' @title Plot Comprehensive GATAI Results
#' @description Create a suite of plots summarizing the effects of GATAI gene removal on phylotranscriptomic patterns.
#'
#' @param phyex_set A PhyloExpressionSet object containing the original gene expression data.
#' @param gatai_result Result list from \code{destroy_pattern()}, containing GATAI analysis output.
#' @param conservation_test Function for conservation test (default: \code{stat_flatline_test}).
#' @param runs_threshold Threshold for gene removal consistency across runs (default: 0.5).
#' @param signature_plot_type Type of signature plot: "separate" for individual plots, "combined" for overlay (default: both options).
#'
#' @return A named list of ggplot/patchwork objects and results:
#'   \item{signature_plots}{Signature plots before/after GATAI and top variance removal}
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
                               conservation_test = stat_flatline_test,
                               runs_threshold = 0.5,
                               signature_plot_type = c("separate", "combined")) {
    check_PhyloExpressionSet(phyex_set)
    plot_type <- match.arg(signature_plot_type)
    removed_genes <- gatai_result$removed_genes   
    n_top_genes <- length(removed_genes)

    # 1.1 Plot TAI signature before/after GATAI
    gatai_set <- remove_genes(phyex_set, removed_genes, new_name = paste(phyex_set@name, "- GATAI removed"))
    q <- 1.0 - n_top_genes/phyex_set@num_genes
    top_var_genes <- genes_top_variance(phyex_set, p = q)
    top_set <- remove_genes(phyex_set, top_var_genes, new_name = paste(phyex_set@name, "- Top variance removed"))

    set.seed(1234)
    all_genes <- phyex_set@gene_ids
    random_genes <- sample(setdiff(all_genes, removed_genes), n_top_genes)
    random_set <- remove_genes(phyex_set, random_genes, new_name = paste(phyex_set@name, "- Random genes removed"))


    if (plot_type == "separate") {
        P1 <- plot_signature(phyex_set, show_p_val = TRUE, conservation_test = conservation_test, colour="blue") +
            ggtitle(paste("Original:", phyex_set@num_genes, "genes"))
        
        P2 <- plot_signature(gatai_set, show_p_val = TRUE, conservation_test = conservation_test, colour="red") +
            ggtitle(paste("GATAI removed:", length(removed_genes), "genes"))
        
        P3 <- plot_signature(top_set, show_p_val = TRUE, conservation_test = conservation_test, colour="purple") +
            ggtitle(paste("Top variance removed:", length(top_var_genes), "genes"))

        P4 <- plot_signature(random_set, show_p_val = TRUE, conservation_test = conservation_test, colour="darkgray") +
            ggtitle(paste("Random genes removed:", length(random_genes), "genes"))
        # Align y-axis scales and combine plots
        y_limits <- range(c(phyex_set@TXI, gatai_set@TXI, top_set@TXI, random_set@TXI)) + c(-0.05, 0.05)
        plots <- list(P1, P2, P3, P4)
        plots <- lapply(plots, \(p) p + ylim(y_limits))
        signature_plots <- (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]])
    }
    else {
        signature_plots <- plot_signature_multiple(c(phyex_set, gatai_set, top_set, random_set), 
                                                   show_p_val = TRUE, 
                                                   conservation_test = conservation_test, 
                                                   colours=c("blue", "red", "purple", "darkgray")) +
            ggtitle(paste(
                "Signature Comparison: Original vs GATAI-removed (", 
                length(removed_genes), "genes) vs Top-variance-removed (", 
                n_top_genes, "genes) vs Random-removed (", 
                length(random_genes), "genes)"
            ))
    }

    # 1.2 Combined signature plot (always include as additional plot)
    signature_combined <- plot_signature_multiple(c(phyex_set, gatai_set), 
                                                 show_p_val = TRUE, 
                                                 conservation_test = conservation_test, 
                                                 colours = c("blue", "red")) +
        ggtitle(paste("Original vs GATAI-removed (", length(removed_genes), "genes)"))


    # 2. Plot removed gene profiles and heatmap
    heatmap_plot <- plot_gene_heatmap(phyex_set, genes=removed_genes, cluster_rows=TRUE, show_gene_ids=TRUE) +
        ggtitle(paste("Expression Heatmap: GATAI-removed genes (", length(removed_genes), "genes)"))
    profiles_plot <- plot_gene_profiles(phyex_set, genes=removed_genes, max_genes=200, transformation="none", colour_by="strata") +
        ggtitle(paste("Gene Expression Profiles: GATAI-removed genes"))
    profiles_plot_facet <- plot_gene_profiles(phyex_set, genes=removed_genes, max_genes=1000, 
                                              facet_by_strata = TRUE, transformation="none", colour_by="strata") +
        ggtitle(paste("Gene Expression Profiles by Phylostrata: GATAI-removed genes")) +
        theme(strip.text = element_text(size = 6))

    gene_space_plot <- plot_gene_space(phyex_set, genes=removed_genes, colour_by="strata") +
        ggtitle(paste("Gene Space Plot: GATAI-removed genes (", length(removed_genes), "genes)"))

    # 3. Mean-variance plot highlighting GATAI-removed genes
    mean_var_plot <- plot_mean_var(phyex_set, highlight_genes = removed_genes) +
        ggtitle(paste("Mean-Variance Plot with GATAI-removed genes highlighted (", length(removed_genes), "genes)"))
    
    # 4. Plot contribution for original and GATAI sets side by side
    contrib_original <- plot_contribution(phyex_set) +
        ggtitle("Gene Contribution: Original Set") +
        theme(legend.position = "none")

    contrib_gatai <- plot_contribution(gatai_set) +
        ggtitle("Gene Contribution: GATAI Set")

    contribution_plots <- contrib_original | contrib_gatai

    # 5. Distribution of all genes by age (raw counts) and GATAI-removed genes (log obs/exp)
    strata_plot_all <- plot_distribution_strata(phyex_set@strata,
                                                as_log_obs_exp = FALSE) +
        ggtitle("Phylostrata Distribution: All Genes")

    strata_plot_removed <- plot_distribution_strata(phyex_set@strata,
                                                    selected_gene_ids = removed_genes,
                                                    as_log_obs_exp = TRUE) +
        ggtitle("Phylostrata Distribution: GATAI-removed genes (Log Obs/Exp)")

    strata_plot <- strata_plot_all / strata_plot_removed
    
    # 6. Compare p values (with individual runs if available)
    original_test <- conservation_test(phyex_set, plot_result = FALSE)
    gatai_test <- conservation_test(gatai_set, plot_result = FALSE)
    
    null_sample <- original_test@null_sample

    # Add top variance and random gene set conservation tests
    top_test <- conservation_test(top_set, plot_result = FALSE)
    random_test <- conservation_test(random_set, plot_result = FALSE)

    test_stats <- data.frame(
        label = c("Original", "GATAI removed", "Top variance removed", "Random genes removed"),
        stat = c(original_test@test_stat, gatai_test@test_stat, top_test@test_stat, random_test@test_stat),
        p_value = c(original_test@p_value, gatai_test@p_value, top_test@p_value, random_test@p_value),
        color = c("blue", "red", "purple", "darkgray")
    )

    # Add individual runs if available
    if ("runs" %in% names(gatai_result)) {
        run_tests <- lapply(gatai_result$runs, function(run_genes) {
            run_set <- remove_genes(phyex_set, run_genes)
            conservation_test(run_set, plot_result = FALSE)
        })
        
        run_stats <- data.frame(
            stat = sapply(run_tests, function(t) t@test_stat),
            p_value = sapply(run_tests, function(t) t@p_value)
        )
        
        # Add run statistics to test_stats for visualization
        run_labels <- paste("Run", seq_len(nrow(run_stats)))
        additional_stats <- data.frame(
            label = run_labels,
            stat = run_stats$stat,
            p_value = run_stats$p_value,
            color = rep("lightgray", nrow(run_stats))
        )
        
        test_stats <- rbind(test_stats, additional_stats)
    }

    null_dist_plot <- ggplot(data.frame(x = null_sample), aes(x = x)) +
        geom_histogram(aes(y = after_stat(density)), bins = 100, fill = "gray67", alpha = 0.7, colour = "gray66") +
        stat_function(fun = original_test@fitting_dist@pdf, args = original_test@params, colour = "gray40") +
        geom_vline(data = test_stats, aes(xintercept = stat, colour = label), linewidth = 1) +
        scale_colour_manual(name = NULL, values = c(
            "Original" = "blue",
            "GATAI removed" = "red", 
            "Top variance removed" = "purple",
            "Random genes removed" = "darkgray",
            setNames(rep("lightgray", sum(grepl("^Run", test_stats$label))), 
                     test_stats$label[grepl("^Run", test_stats$label)])
        )) +
        labs(x = "Score", y = "Density",
             title = "Conservation Test Comparison",
             subtitle = paste("Test:", original_test@method_name)) +
        annotate(
            "text",
            x = test_stats$stat - 0.05 * diff(range(null_sample)),
            y = max(density(null_sample)$y) * seq(0.9, 0.1, length.out = nrow(test_stats)),
            label = sapply(test_stats$p_value, function(p) formatC(p, format = "e", digits = 2)),
            colour = test_stats$color,
            hjust = 1,
            size = 3.5
        ) +
        theme_minimal()

    result_list <- list(
        signature_plots = signature_plots,
        signature_combined = signature_combined,
        mean_var_plot = mean_var_plot,
        contribution_plots = contribution_plots,
        strata_plot = strata_plot,
        heatmap_plot = heatmap_plot,
        profiles_plot = profiles_plot,
        profiles_plot_facet = profiles_plot_facet,
        gene_space_plot = gene_space_plot,
        null_dist_plot = null_dist_plot
    )

    # 7. (optional) Convergence plots

    # check if runs is included in the result
    if ("runs" %in% names(gatai_result)) {
        message("Computing GATAI result convergence plots. This may take a while...")
        result_list$convergence_plots <- full_gatai_convergence_plot(phyex_set, gatai_result$runs, p=runs_threshold) + 
            theme_minimal(base_size = 2) + 
            theme(plot.title = element_text(size = 6),
                  strip.text = element_text(size = 4))
        
        # 8. Individual runs signature comparison plot
        message("Creating individual runs signature comparison plot...")
        
        # Create sets for each individual run
        run_sets <- lapply(seq_along(gatai_result$runs), function(i) {
            run_genes <- gatai_result$runs[[i]]
            remove_genes(phyex_set, run_genes, new_name = paste("Run", i, "removed"))
        })
        
        # Combine original, GATAI consensus, and individual runs
        all_sets <- c(list(phyex_set, gatai_set), run_sets)
        
        # Create colour palette: blue for original, red for GATAI, light colours for runs
        n_runs <- length(gatai_result$runs)
        run_colours <- rep("lightgray", n_runs)
        all_colours <- c("blue", "red", run_colours)
        
        result_list$individual_runs_signature <- plot_signature_multiple(
            all_sets,
            show_p_val = FALSE,  # Too many p-values would clutter the plot
            conservation_test = conservation_test,
            colours = all_colours
        ) +
            ggtitle(paste("Signature Comparison: Original vs GATAI Consensus vs Individual Runs (", 
                         n_runs, "runs)")) +
            theme(legend.position = "bottom") +
            guides(colour = guide_legend(nrow = 2))
    }

    
    return(result_list)
}

#' @title Save GATAI Analysis Results to PDF
#' @description Save removed gene IDs and all GATAI analysis plots to a PDF file.
#'
#' @param phyex_set A PhyloExpressionSet object containing the original gene expression data.
#' @param gatai_result Result list from \code{destroy_pattern()}, containing GATAI analysis output.
#' @param analysis_dir Directory to save the PDF file.
#' @param prefix Optional prefix for the PDF filename (default: "report").
#' @param ... Additional arguments passed to \code{plot_gatai_results()}.
#'
#' @return Invisibly returns the path to the saved PDF.
#'
#' @examples
#' # Save results after running destroy_pattern
#' # save_gatai_results_pdf("results/", phyex_set, gatai_result)
#'
#' @importFrom stats density
#' @importFrom graphics plot.new title text
#' @export
save_gatai_results_pdf <- function(phyex_set,
                                   gatai_result,
                                   analysis_dir = "gatai_analysis",
                                   prefix = "report",
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