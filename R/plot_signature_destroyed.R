#' @title Gene removal analysis using GATAI and comparison with random removal
#' @description Function to analyze the effect of GATAI gene removal on phylotranscriptomic signatures compared to random gene removal.
#' @param phyex_set a PhyloExpressionSet object.
#' @param plot_type visualization type: "separate" for individual plots or "combined" for overlay plot.
#' @param conservation_test the conservation test to use for p-value calculation.
#' @param n_random_removal number of randomly removed genes. Default is "gatai" (same number as GATAI removed).
#' @param ... additional arguments passed to destroy_pattern function.
#' @author Filipa Martins Costa
#' @export
plot_signature_destroyed <- function(phyex_set,
                                      plot_type = "separate",
                                      conservation_test = flatline_test,
                                      n_random_removal = "gatai",
                                      ...) {
    if (!plot_type %in% c("separate", "combined")) {
        stop("Plot type '", plot_type, "' is not available for this function. Please specify a plot type supported by this function.", call. = FALSE)
    }

    # Get GATAI removed genes
    gatai_result <- destroy_pattern(phyex_set, ...)
    removed_genes <- gatai_result$genes_removed

    # Determine number of random genes to remove
    if (is.character(n_random_removal) && n_random_removal == "gatai") {
        n_random_genes <- length(removed_genes)
    } else if (is.numeric(n_random_removal)) {
        n_random_genes <- n_random_removal
    } else {
        message("Invalid value for 'n_random_removal'. It should be 'gatai' or a numeric value.")
        message("Removing the same number as GATAI.")
        n_random_genes <- length(removed_genes)
    }

    # Create the three datasets
    original_set <- phyex_set
    gatai_set <- remove_genes(phyex_set, removed_genes, new_name = paste(phyex_set@name, "- GATAI removed"))

    random_removed_genes <- sample(phyex_set@gene_ids, n_random_genes)
    random_set <- remove_genes(phyex_set, random_removed_genes, new_name = paste(phyex_set@name, "- Random removed"))

    if (plot_type == "separate") {
        # Create individual plots
        P1 <- plot_signature(original_set, show_p_val = TRUE, conservation_test = conservation_test) +
            ggtitle(paste("Original:", original_set@num_genes, "genes"))

        P2 <- plot_signature(gatai_set, show_p_val = TRUE, conservation_test = conservation_test) +
            ggtitle(paste("GATAI removed:", length(removed_genes), "genes"))

        P3 <- plot_signature(random_set, show_p_val = TRUE, conservation_test = conservation_test) +
            ggtitle(paste("Random removed:", n_random_genes, "genes"))

        # Align y-axis scales
        y_data <- c(original_set@TXI, gatai_set@TXI, random_set@TXI)
        y_limits <- c(min(y_data) - 0.05, max(y_data) + 0.05)

        P1 <- P1 + ylim(y_limits)
        P2 <- P2 + ylim(y_limits)
        P3 <- P3 + ylim(y_limits)

        plots <- cowplot::plot_grid(P1, P2, P3, ncol = 3)
        print(plots)
    }

    if (plot_type == "combined") {
        phyex_sets <- list(original_set, gatai_set, random_set)

        p <- plot_signature_multiple(phyex_sets,
            legend_title = "Dataset",
            show_p_val = TRUE,
            conservation_test = conservation_test
        )

        print(p)
    }

    return(list(
        gatai_removed_genes = removed_genes,
        random_removed_genes = random_removed_genes,
        original_set = original_set,
        gatai_set = gatai_set,
        random_set = random_set
    ))
}
