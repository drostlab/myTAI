#' @title Partial TAI Distribution Plotting Functions
#' @description Functions for plotting and comparing partial TAI distributions using PhyloExpressionSet S7 objects.
#' @import ggplot2
#' @importFrom dplyr inner_join mutate filter bind_rows left_join
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom ggridges geom_density_ridges
#' @importFrom cowplot plot_grid
#' @importFrom purrr map
#' @importFrom stats ks.test qqplot

#' @title Comparing expression/partial TAI distributions across developmental stages
#' @description \emph{plot_distribution_pTAI} generates 2 plots that help to compare the distribution
#' of the quotient of expression by partial TAI through various developmental stages or cell types, highlighting each stage with 
#' distinct colors.
#' @param phyex_set A PhyloExpressionSet object (BulkPhyloExpressionSet or ScPhyloExpressionSet).
#' @param stages A numeric vector specifying the indices of the stages to compare. Each index 
#' corresponds to a stage in the PhyloExpressionSet. If NULL, all stages are used.
#' @param xlab Label of x-axis.
#' @param ylab Label of y-axis.
#' @param main Figure title.
#' @author Filipa Martins Costa
#' @export
plot_distribution_pTAI <- function(phyex_set,
                                   stages = NULL,
                                   xlab = "Expression / Partial TAI",
                                   ylab = "Density",
                                   main = "Density Distribution of Expression / Partial TAI by Developmental Stage") {
    
    # Use all stages if none specified
    if (is.null(stages)) {
        stages <- 1:phyex_set@num_identities
    }
    
    if (any(stages > phyex_set@num_identities)) {
        stop("Some indices in 'stages' exceed the number of identities in the PhyloExpressionSet.")
    }
    
    # Get partial TAI matrix
    partial_TAI_matrix <- pTXI(phyex_set)[, stages, drop = FALSE]
    
    # Get expression data for selected stages
    expression_matrix <- phyex_set@expression_collapsed[, stages, drop = FALSE]
    
    # Convert to long format
    partial_TAI_df <- tibble::rownames_to_column(as.data.frame(partial_TAI_matrix), var = "GeneID")
    partial_TAI_long <- tidyr::pivot_longer(
        partial_TAI_df,
        cols = -GeneID,
        names_to = "Stage",
        values_to = "PartialTAI"
    )
    
    expression_df <- tibble::rownames_to_column(as.data.frame(expression_matrix), var = "GeneID")
    expression_long <- tidyr::pivot_longer(
        expression_df,
        cols = -GeneID,
        names_to = "Stage",
        values_to = "Expression"
    )
    
    # Combine data
    combined_long <- dplyr::inner_join(partial_TAI_long, expression_long, by = c("GeneID", "Stage"))
    combined_long <- dplyr::mutate(combined_long, Ratio = Expression / PartialTAI)
    
    # Set factor levels for consistent ordering
    combined_long$Stage <- factor(combined_long$Stage, levels = colnames(partial_TAI_matrix))
    
    # Generate colors
    qual_col_pals <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    col_vector <- unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    set.seed(123)
    colors <- sample(col_vector, ncol(partial_TAI_matrix))
    
    # Create density plot
    P1 <- ggplot2::ggplot(combined_long, ggplot2::aes(x = Ratio, fill = Stage)) +
        ggplot2::geom_density(alpha = 0.7, color = "black") +
        ggplot2::labs(
            x = xlab,
            y = ylab,
            title = phyex_set@name
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
            axis.title = ggplot2::element_text(size = 12)
        ) +
        ggplot2::scale_fill_manual(values = colors)
    
    # Create ridgeline plot
    P2 <- ggplot2::ggplot(combined_long, ggplot2::aes(x = Ratio, y = Stage, fill = Stage)) +
        ggridges::geom_density_ridges(alpha = 0.7, color = "black", scale = 1) + 
        ggplot2::labs(
            x = xlab,
            y = phyex_set@identities_label,
            title = phyex_set@name
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
            axis.title = ggplot2::element_text(size = 12)
        ) +
        ggplot2::scale_fill_manual(values = colors)
    
    cowplot::plot_grid(P1, P2, labels = main)
}

#' @title QQ plot comparing partial TAI distributions across developmental stages against a reference stage
#' @description \emph{plot_distribution_partialTAI_qqplot} generates a QQ plot to compare the partial TAI 
#' distributions of various developmental stages against a reference stage.
#' It visualizes quantile differences between the reference and other stages,
#' highlights each stage with distinct colors, and annotates the plot with the p-values
#' from the nonparametric \code{\link{ks.test}} to indicate the significance of distribution differences.
#' @param phyex_set A PhyloExpressionSet object (BulkPhyloExpressionSet or ScPhyloExpressionSet).
#' @param reference_stage_index An integer specifying the index of the reference developmental stage. 
#' The partial TAI distribution of this stage will be used as the reference 
#' for comparisons with other stages (default: stage index 1).
#' @param xlab Label of x-axis.
#' @param ylab Label of y-axis.
#' @param main Figure title.
#' @param alpha Transparency of the points.
#' @param size Size of the points.
#' @author Filipa Martins Costa
#' @export
plot_distribution_pTAI_qqplot <- function(phyex_set, 
                                               reference_stage_index = 1,
                                               xlab = "Quantiles of Reference Stage",
                                               ylab = "Quantiles of Other Stages",
                                               main = "QQ Plot: Developmental Stages vs Reference Stage (p-values from KS test)",
                                               alpha = 0.7,
                                               size = 1.2) {
    
    partial_TAI_matrix <- pTXI(phyex_set)
    
        if (reference_stage_index > phyex_set@num_identities) {
        stop("The specified reference stage index exceeds the number of identities in the PhyloExpressionSet. Please provide a valid index.")
    }
    
    # Define the reference stage
    reference_stage <- partial_TAI_matrix[, reference_stage_index]
    
    # Function to compute QQ data and KS p-value for a given stage
    compute_qq_ks <- function(stage_index) {
        stage_values <- partial_TAI_matrix[, stage_index]
        
        # Compute QQ plot data
        qq <- qqplot(reference_stage, stage_values, plot.it = FALSE)
        
        # Perform KS test
        ks_p_value <- ks.test(reference_stage, stage_values)$p.value
        
        # Return a list containing QQ data and KS p-value
        list(
            qq_data = data.frame(
                reference = qq$x,  # quantiles of the reference stage
                stage_quantiles = qq$y,  # quantiles of the current stage
                stage = as.character(phyex_set@identities)[stage_index]  # stage name
            ),
            ks_result = data.frame(
                stage = as.character(phyex_set@identities)[stage_index],
                p_value = ks_p_value
            )
        )
    }
    
    # Apply the function to all stages except the reference
        results <- purrr::map(setdiff(1:phyex_set@num_identities, reference_stage_index), compute_qq_ks)
    
    # Combine QQ plot data
    qq_data <- dplyr::bind_rows(purrr::map(results, "qq_data"))
    
    # Combine KS test results
    ks_results <- dplyr::bind_rows(purrr::map(results, "ks_result"))
    
    # Add p-values to stage labels
    ks_results <- dplyr::mutate(
        ks_results,
        stage_label = paste(stage, "(p-value =", format(p_value, digits = 3), ")")
    )
    
    # Merge stage labels into the qq_data frame
    qq_data <- dplyr::left_join(
        qq_data,
        ks_results[, c("stage", "stage_label")],
        by = "stage"
    )
    qq_data <- dplyr::mutate(
        qq_data,
        stage_label = factor(stage_label, levels = unique(ks_results$stage_label))
    )
    
    ggplot2::ggplot(qq_data, ggplot2::aes(x = reference, y = stage_quantiles, color = stage_label)) +
        ggplot2::geom_point(alpha = alpha, size = size) +  
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +  
        ggplot2::labs(
            title = paste(phyex_set@name, "-", main),
            x = xlab,
            y = ylab,
            color = paste(phyex_set@identities_label, "(p-value)")
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "top") +
        ggplot2::scale_color_brewer(palette = "Set1")
}
