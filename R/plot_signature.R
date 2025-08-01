#' @title Plot Transcriptomic Signature
#' @description Create a plot of the transcriptomic index signature across developmental stages
#' or cell types, with options for showing individual samples/cells and statistical testing.
#' 
#' @param phyex_set A PhyloExpressionSet object (BulkPhyloExpressionSet or ScPhyloExpressionSet)
#' @param show_reps Logical, whether to show individual replicates
#' @param show_p_val Logical, whether to show conservation test p-value
#' @param conservation_test Function, conservation test to use for p-value calculation
#' @param colour Character, custom color for the plot elements
#' @param ... Additional arguments passed to specific methods
#' 
#' @return A ggplot2 object showing the transcriptomic signature
#' 
#' @details
#' This function creates visualizations appropriate for the data type:
#' 
#' **Bulk data (BulkPhyloExpressionSet):**
#' - Line plots showing TXI trends across developmental stages
#' - Optional individual biological replicates as jittered points
#' - Optional conservation test p-values
#' 
#' **Single-cell data (ScPhyloExpressionSet):**
#' - Violin plots showing TXI distributions across cell types
#' - Mean TXI values overlaid as line
#' - Optional individual cells using geom_sina for better visualization
#' 
#' @examples
#' # Basic signature plot for bulk data
#' # p1 <- plot_signature(bulk_phyex_set)
#' 
#' # Bulk plot with replicates and p-value
#' # p2 <- plot_signature(bulk_phyex_set, show_reps = TRUE, show_p_val = TRUE)
#' 
#' # Single-cell plot with individual cells
#' # p3 <- plot_signature(sc_phyex_set, show_reps = TRUE)
#' 
#' @import ggplot2
#' @export
plot_signature <- S7::new_generic("plot_signature", "phyex_set",
    function(phyex_set,
             show_reps = TRUE,
             show_p_val = TRUE,
             conservation_test = stat_flatline_test,
             colour = NULL,
             ...) {
        S7::S7_dispatch()
    }
)

#' @export
S7::method(plot_signature, BulkPhyloExpressionSet) <- function(phyex_set,
                                                               show_reps = TRUE,
                                                               show_p_val = TRUE,
                                                               conservation_test = stat_flatline_test,
                                                               colour = NULL,
                                                               show_bootstraps = FALSE,
                                                               show_stddev = TRUE,
                                                               ...) {
                                                               
                                                               args <- list(...)
    
    # Create main TXI line
    df_main <- tibble::tibble(
        Identity = phyex_set@identities,
        TXI = phyex_set@TXI
    )
    
    # Start with line plot
    p <- ggplot()

    # Optionally plot bootstrapped TXI lines
    if (show_bootstraps) {
        df_b <- as_tibble(phyex_set@bootstrapped_txis) |>
            rowid_to_column("Id") |>
            tidyr::pivot_longer(cols=phyex_set@identities, names_to="Identity", values_to = "TXI")
        p <- p + geom_line(data=df_b,
                           aes(x=Identity, y=TXI, group=Id, colour=phyex_set@name),
                           alpha=0.01)
    }

    p <- p +
        geom_line(data=df_main, aes(x = Identity, y = TXI, group = 1),
                  colour = "black", lwd = 2.3, lineend = "round") +
        geom_line(data=df_main, 
                  aes(x = Identity, y = TXI, group = 1, colour = phyex_set@name), 
                  lwd = 1.5, lineend = "round")


    # Add replicate dots and CI ribbon if requested
    if (show_reps) {
        df_samples <- tibble::tibble(
            Identity = phyex_set@groups,
            TXI = phyex_set@TXI_sample
        )

        # Compute mean and SD for each identity using bootstrapped TXI values
        std_dev <- TXI_std_dev(phyex_set)

        # Build a tibble for ribbon: mean ± 1 SD
        df_sd <- tibble::tibble(
            Identity = factor(phyex_set@identities, levels = unique(as.character(phyex_set@identities))),
            lb = phyex_set@TXI - std_dev,
            ub = phyex_set@TXI + std_dev
        )

        if (show_stddev)
            p <- p +
                geom_ribbon(
                    data = df_sd,
                    aes(x = Identity, ymin = lb, ymax = ub, fill = phyex_set@name, group=1), alpha = 0.3, inherit.aes = FALSE
                )
            
        p <- p +
            geom_jitter(
                data = df_samples,
                aes(x = Identity, y = TXI, fill = phyex_set@name),
                shape = 21, colour = "black", size = 1.5, stroke = 0.5, 
                width = 0.05, alpha = 0.7
            )
    }

    # Add labels and theme
    p <- p +
        labs(
            x = phyex_set@identities_label,
            y = phyex_set@index_full_name
        ) +
        guides(colour = "none", fill = "none") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

    # Show p value for conservation tests
    if (show_p_val) {
        if ("modules" %in% names(args))
            t <- conservation_test(phyex_set, plot_result = FALSE, modules=args$modules)
        else
            t <- conservation_test(phyex_set, plot_result = FALSE)
        label <- exp_p(t@p_value)
        p <- p +
            annotate("text",
                label = label, parse = TRUE,
                x = phyex_set@num_identities * 0.7, 
                y = mean(phyex_set@TXI_sample) + 0.1,
                size = 6 # Increased text size
            )
    }

    # Apply custom colour if specified
    if (!is.null(colour)) {
        p <- p +
            scale_color_manual(values = c(colour)) +
            scale_fill_manual(values = c(colour))
    }

    return(p)
}

#' @export
S7::method(plot_signature, ScPhyloExpressionSet) <- function(phyex_set, 
                                                             show_reps = TRUE,
                                                             colour = NULL,
                                                             show_p_val = TRUE,
                                                             conservation_test = stat_flatline_test,
                                                             ...) {
    
    # Prepare data for plotting
    df_samples <- tibble::tibble(
        Identity = phyex_set@groups,
        TXI = phyex_set@TXI_sample
    )
    
    # Use the proper TXI values (aggregated values per identity)
    df_main <- tibble::tibble(
        Identity = phyex_set@identities,
        TXI = phyex_set@TXI
    )
    
    # Create base plot
    if (show_reps) {
        # Show violin plots with individual cells as reps
        p <- ggplot(df_samples, aes(x = Identity, y = TXI, fill = Identity)) +
            geom_violin(alpha = 0.7, scale = "width") +
            ggforce::geom_sina(
                size = 0.1, alpha = 0.5,
                colour = "black"
            )
        
        # Add TXI line plot on top
        p <- p + geom_line(
            data = df_main,
            aes(x = Identity, y = TXI, group = 1),
            colour = "black", lwd = 2.3, lineend = "round", inherit.aes = FALSE
        ) +
            geom_line(
                data = df_main,
                aes(x = Identity, y = TXI, group = 1, colour = phyex_set@name),
                lwd = 1.5, lineend = "round", inherit.aes = FALSE
            )
        
        # Add mean points (proper TXI values) on top
        p <- p + geom_point(
            data = df_main,
            aes(x = Identity, y = TXI),
            shape = 21, colour = "black", size = 2, stroke = 0.8,
            fill = "white", inherit.aes = FALSE
        )
    } else {
        # Simple plot showing only aggregated TXI values with line
        p <- ggplot(df_main, aes(x = Identity, y = TXI, fill = Identity, group = 1)) +
            geom_line(colour = "black", lwd = 2.3, lineend = "round") +
            geom_line(aes(colour = phyex_set@name), lwd = 1.5, lineend = "round") +
            geom_point(size = 3, shape = 21, colour = "black", stroke = 0.8)
    }
    
    p <- p +
        labs(
            x = phyex_set@identities_label,
            y = phyex_set@index_full_name
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none"
        )
    
    # Show p value for conservation tests
    if (show_p_val) {
        t <- conservation_test(phyex_set, plot_result = FALSE)
        label <- exp_p(t@p_value)
        p <- p +
            annotate("text",
                label = label, parse = TRUE,
                x = phyex_set@num_identities * 0.7, 
                y = max(df_samples$TXI, na.rm = TRUE) * 0.9
            )
    }
    
    # Apply colour scheme
    if (!is.null(colour)) {
        p <- p +
            scale_fill_manual(values = rep(colour, phyex_set@num_identities))
    } else {
        p <- p + scale_fill_viridis_d()
    }
    
    return(p)
}