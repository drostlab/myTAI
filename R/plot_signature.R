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
#' - Sina plots showing TXI distributions across cell types or other identities
#' - Mean TXI values overlaid as line
#' - Optional individual cells using geom_sina for better visualization
#' - Flexible identity selection from metadata via additional parameters:
#'   - `primary_identity`: Character, name of metadata column for x-axis (default: current selected identities)
#'   - `secondary_identity`: Character, name of metadata column for coloring/faceting
#'   - `facet_by_secondary`: Logical, whether to facet by secondary identity (default: FALSE uses colouring)
#' 
#' @examples
#' phyex_set <- example_phyex_set |>
#'     select_genes(example_phyex_set@gene_ids[1:100])
#' phyex_set@null_conservation_sample_size <- 500
#' 
#' # Basic signature plot for bulk data
#' p1 <- plot_signature(example_phyex_set)
#' 
#' # Bulk plot with replicates and p-value
#' p2 <- plot_signature(example_phyex_set, show_reps = TRUE, show_p_val = TRUE)
#' 
#' phyex_set_sc <- example_phyex_set_sc
#' phyex_set_sc@null_conservation_sample_size <- 500
#' 
#' # Single-cell plot with individual cells
#' p3 <- plot_signature(phyex_set_sc, show_reps = TRUE)
#' 
#' # Single-cell plot with custom primary identity
#' p4 <- plot_signature(phyex_set_sc, primary_identity = "day")
#' 
#' # Single-cell plot with primary and secondary identities (colored)
#' p5 <- plot_signature(phyex_set_sc, 
#'                      primary_identity = "day", 
#'                      secondary_identity = "condition")
#' 
#' # Single-cell plot with faceting by secondary identity
#' p6 <- plot_signature(phyex_set_sc, 
#'                      primary_identity = "day", 
#'                      secondary_identity = "condition", 
#'                      facet_by_secondary = TRUE)

#' @import ggplot2
#' @importFrom dplyr group_by summarise mutate
#' @importFrom tibble tibble
#' @importFrom ggforce geom_sina
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
        p <- p + ggforce::geom_sina(data=df_b,
                           aes(x=Identity, y=TXI, colour=phyex_set@name),
                           alpha=0.3)
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
        if ("modules" %in% names(args)) {
            t <- conservation_test(phyex_set, plot_result = FALSE, modules=args$modules)
        }
        else
            t <- conservation_test(phyex_set, plot_result = FALSE)
        message(paste("Ran", t@method_name))
        if ("modules" %in% names(args)) {
            modules <- args$modules
            message("Modules: \n early = {",paste0(phyex_set@sample_names[modules[[1]]], " "),"}","\n","mid = {",paste0(phyex_set@sample_names[modules[[2]]], " "),"}","\n","late = {",paste0(phyex_set@sample_names[modules[[3]]], " "),"}")
        }
        message("Significance status of signature: ", 
            ifelse(as.numeric(t@p_value) <= 0.05, "significant.", "not significant (= no evolutionary signature in the transcriptome)."))
        label <- exp_p(t@p_value)
        p <- p +
            annotate("text",
                label = label, parse = TRUE,
                x = 1, 
                y = Inf,
                hjust = 0,
                vjust = 1,
                size = 6 
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
                                                             primary_identity = NULL,
                                                             secondary_identity = NULL,
                                                             facet_by_secondary = FALSE,
                                                             ...) {
    
    # Prepare data using helper function
    plot_data <- .prepare_sc_plot_data(phyex_set, primary_identity, secondary_identity)
    df_samples <- plot_data$samples
    df_main <- plot_data$main
    primary_label <- plot_data$primary_label
    secondary_label <- plot_data$secondary_label
    
    # Determine which identity to use for coloring
    color_by_identity <- if (!is.null(secondary_identity)) secondary_identity else primary_identity
    if (is.null(color_by_identity)) color_by_identity <- phyex_set@identities_label
    
    # Check for custom colors and show message if using defaults
    has_custom_colors <- !is.null(phyex_set@idents_colours) && 
                         !is.null(phyex_set@idents_colours[[color_by_identity]])
    if (!has_custom_colors && is.null(colour)) {
        message(sprintf(
            "Using default colors for identity '%s'. To set custom colors, assign to phyex_set@idents_colours[['%s']] <- c(value1 = 'color1', value2 = 'color2', ...)",
            color_by_identity, color_by_identity
        ))
    }
    
    # Handle mixed x-axis for secondary identity without faceting
    if (!is.null(secondary_identity) && !facet_by_secondary) {
        # Create compound x-axis labels: Primary_Secondary
        df_samples <- df_samples |>
            dplyr::mutate(
                Mixed_X = paste(Primary, Secondary, sep = "_"),
                Mixed_X = factor(Mixed_X, levels = unique(Mixed_X[order(Primary, Secondary)]))
            )
        
        # Update main data for line plot (average within each Mixed_X group)
        df_main <- df_samples |>
            dplyr::group_by(Mixed_X, Primary, Secondary) |>
            dplyr::summarise(TXI = mean(TXI, na.rm = TRUE), .groups = "drop")
    }
    
    # Create base plot
    if (show_reps) {
        if (!is.null(secondary_identity) && !facet_by_secondary) {
            # Mixed x-axis with coloring by secondary identity
            p <- ggplot(df_samples, aes(x = Mixed_X, y = TXI, color = Secondary)) +
                # Add individual points with sina
                ggforce::geom_sina(size = 0.8, alpha = 0.7) +
                # Add faint black line connecting means
                stat_summary(fun = mean, geom = "line", 
                           color = "black", alpha = 0.4, size = 0.8, 
                           aes(group = 1), show.legend = FALSE) +
                # Add prominent mean lines for each identity with matching colors
                geom_boxplot(width=0.5, outlier.shape=NA, color="black",  fill="white", alpha=0.5)
        } else {
            # Standard plot - color by primary identity
            p <- ggplot(df_samples, aes(x = Primary, y = TXI, color = Primary)) +
                # Add individual points with sina
                ggforce::geom_sina(size = 0.8, alpha = 0.7) +
                # Add faint black line connecting means
                stat_summary(fun = mean, geom = "line", 
                           color = "black", alpha = 0.4, size = 0.8, 
                           aes(group = 1), show.legend = FALSE) +
                # Add prominent mean lines for each identity with matching colors
                geom_boxplot(width=0.5, outlier.shape=NA, color="black",  fill="white", alpha=0.5)
        }
    } else {
        # Simple plot showing only aggregated TXI values with line
        if (!is.null(secondary_identity) && !facet_by_secondary) {
            p <- ggplot(df_main, aes(x = Mixed_X, y = TXI, group = 1)) +
                geom_line(colour = "black", lwd = 2.3, lineend = "round") +
                geom_line(aes(colour = phyex_set@name), lwd = 1.5, lineend = "round") +
                geom_point(size = 3, shape = 21, colour = "black", stroke = 0.8, fill = "white")
        } else {
            p <- ggplot(df_main, aes(x = Primary, y = TXI, group = 1)) +
                geom_line(colour = "black", lwd = 2.3, lineend = "round") +
                geom_line(aes(colour = phyex_set@name), lwd = 1.5, lineend = "round") +
                geom_point(size = 3, shape = 21, colour = "black", stroke = 0.8, fill = "white")
        }
    }
    
    # Add faceting if requested
    if (!is.null(secondary_identity) && facet_by_secondary) {
        p <- p + facet_wrap(~ Secondary, scales = "free_y")
    }
    
    # Add labels and theme
    x_label <- if (!is.null(secondary_identity) && !facet_by_secondary) {
        paste(primary_label, secondary_label, sep = " / ")
    } else {
        primary_label
    }
    
    p <- p +
        labs(
            x = x_label,
            y = phyex_set@index_full_name
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = if(!is.null(secondary_identity) && !facet_by_secondary) "right" else "none"
        )
    
    # Hide legend for the dataset name when show_reps = FALSE (consistent with bulk)
    if (!show_reps) {
        p <- p + guides(colour = "none", fill = "none")
    }
    
    # Show p value for conservation tests (only if using default identities)
    if (show_p_val && is.null(primary_identity)) {
        t <- conservation_test(phyex_set, plot_result = FALSE)
        label <- exp_p(t@p_value)
        p <- p +
            annotate("text",
                label = label, parse = TRUE,
                x = 1, 
                y = Inf,
                hjust = 0, vjust = 1,
                size = 6
            )
    }
    
    # Apply colour scheme
    if (!is.null(colour)) {
        # Manual color override
        if (!is.null(secondary_identity) && !facet_by_secondary) {
            n_secondary <- length(unique(df_samples$Secondary))
            p <- p + scale_color_manual(values = rep(colour, n_secondary))
        } else {
            n_primary <- length(unique(df_samples$Primary))
            p <- p + scale_color_manual(values = rep(colour, n_primary))
        }
    } else if (has_custom_colors) {
        # Use custom colors from the object - ensure they match the data levels
        custom_colors <- phyex_set@idents_colours[[color_by_identity]]
        
        if (!is.null(secondary_identity) && !facet_by_secondary) {
            actual_levels <- levels(factor(df_samples$Secondary))
            if (all(actual_levels %in% names(custom_colors))) {
                p <- p + scale_color_manual(values = custom_colors[actual_levels], name = color_by_identity)
            } else {
                # Fallback to viridis if mismatch
                p <- p + scale_color_viridis_d(name = color_by_identity)
            }
        } else {
            actual_levels <- levels(factor(df_samples$Primary))
            if (all(actual_levels %in% names(custom_colors))) {
                p <- p + scale_color_manual(values = custom_colors[actual_levels], name = color_by_identity)
            } else {
                # Fallback to viridis if mismatch
                p <- p + scale_color_viridis_d(name = color_by_identity)
            }
        }
    } else {
        # Use default viridis colors
        p <- p + scale_color_viridis_d(name = color_by_identity)
    }
    
    return(p)
}

## HELPER FUNCTIONS

#' @title Prepare Single-Cell Plot Data
#' @description Prepare data frame for single-cell plotting with flexible identities
#' @param phyex_set A ScPhyloExpressionSet object
#' @param primary_identity Primary identity column name
#' @param secondary_identity Secondary identity column name (optional)
#' @return List with sample data and aggregated data
#' @keywords internal
.prepare_sc_plot_data <- function(phyex_set, primary_identity = NULL, secondary_identity = NULL) {
    # Get metadata
    metadata <- phyex_set@metadata
    
    # Determine primary identity
    if (is.null(primary_identity)) {
        # Use current selected identities via groups
        primary_values <- phyex_set@groups
        primary_label <- phyex_set@identities_label
    } else {
        # Use specified column from metadata
        if (!primary_identity %in% colnames(metadata)) {
            stop(sprintf("Primary identity '%s' not found in metadata. Available columns: %s",
                        primary_identity, paste(colnames(metadata), collapse = ", ")))
        }
        primary_values <- metadata[[primary_identity]]
        primary_label <- primary_identity
    }
    
    # Create base data frame with TXI values
    df_samples <- tibble::tibble(
        Primary = primary_values,
        TXI = phyex_set@TXI_sample
    )
    
    # Add secondary identity if specified
    if (!is.null(secondary_identity)) {
        if (!secondary_identity %in% colnames(metadata)) {
            stop(sprintf("Secondary identity '%s' not found in metadata. Available columns: %s",
                        secondary_identity, paste(colnames(metadata), collapse = ", ")))
        }
        df_samples$Secondary <- metadata[[secondary_identity]]
        secondary_label <- secondary_identity
    } else {
        secondary_label <- NULL
    }
    
    # Calculate aggregated TXI values for line plot
    if (is.null(secondary_identity)) {
        # Simple aggregation by primary identity only
        df_main <- df_samples |>
            dplyr::group_by(Primary) |>
            dplyr::summarise(TXI = mean(TXI, na.rm = TRUE), .groups = "drop")
    } else {
        # Aggregate by both primary and secondary, then average over secondary for line
        df_main <- df_samples |>
            dplyr::group_by(Primary, Secondary) |>
            dplyr::summarise(TXI = mean(TXI, na.rm = TRUE), .groups = "drop") |>
            dplyr::group_by(Primary) |>
            dplyr::summarise(TXI = mean(TXI, na.rm = TRUE), .groups = "drop")
    }
    
    return(list(
        samples = df_samples,
        main = df_main,
        primary_label = primary_label,
        secondary_label = secondary_label
    ))
}