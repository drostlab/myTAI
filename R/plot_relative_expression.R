#' @title Relative Expression Functions
#' @description Functions for computing and plotting relative expression profiles using PhyloExpressionSet S7 objects.
#' @import ggplot2
#' @importFrom dplyr filter case_when
#' @importFrom reshape2 melt
#' @importFrom stats kruskal.test
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales pretty_breaks

# Helper Functions --------------------------------------------------------

#' @title Transform to Relative Expression Levels
#' @description Computes the relative expression profile for a given gene expression matrix. 
#' The relative expression is calculated by normalizing the column means of the matrix to a [0, 1] scale.
#' @param count_matrix A numeric matrix where columns represent developmental stages and rows represent genes.
#' @return A numeric vector of relative expression values for each stage (column) in the input matrix.
#' @export
relative_expression <- function(count_matrix) {
    col_means <- colMeans(count_matrix)
    f_max <- max(col_means)
    f_min <- min(col_means)
    if (f_max == f_min) return(rep(0, length(col_means)))
    return((col_means - f_min) / (f_max - f_min))
}

#' @title Compute Relative Expression Matrix for PhyloExpressionSet
#' @description Computes relative expression profiles for all age categories in a PhyloExpressionSet.
#' @param phyex_set A PhyloExpressionSet object.
#' @return A matrix with age categories as rows and conditions as columns, containing relative expression values.
#' @export
rel_exp_matrix <- function(phyex_set) {
    if (phyex_set@num_conditions < 2) stop("You need at least 2 conditions to compute relative expression levels.")
    
    age_vec <- as.integer(phyex_set@strata)
    counts <- phyex_set@counts_collapsed
    
    # Compute relative expression for each age category
    rel_exp_mat <- t(sapply(sort(unique(age_vec)), function(a) {
        idx <- which(age_vec == a)
        mat <- counts[idx, , drop = FALSE]
        relative_expression(mat)
    }))
    
    rownames(rel_exp_mat) <- sort(unique(age_vec))
    colnames(rel_exp_mat) <- as.character(phyex_set@conditions)
    return(rel_exp_mat)
}

# Plotting Functions ------------------------------------------------------

#' @title Plot Relative Expression Profiles (Line Plot)
#' @description Plots relative expression profiles for age categories using a PhyloExpressionSet S7 object.
#' @param phyex_set A PhyloExpressionSet object.
#' @param groups A list of integer vectors specifying age categories (e.g., phylostrata) for each group (1 or 2 groups).
#' @param modules Optional list for shading modules: list(early=..., mid=..., late=...).
#' @param adjust_range Logical, adjust y-axis range for both panels (if 2 groups).
#' @param alpha Transparency for shaded module area.
#' @param ... Further arguments passed to ggplot2 geoms.
#' @return ggplot2 object or list of ggplot2 objects.
#' @export
plot_relative_expression_line <- function(
    phyex_set,
    groups,
    modules = NULL,
    adjust_range = TRUE,
    alpha = 0.1,
    ...) {
    if (is.null(groups) || !is.list(groups) || length(groups) < 1) stop("groups must be a non-empty list.")
    if (length(groups) > 2) stop("For line plot, specify at most 2 groups.")

    age_vec <- as.integer(phyex_set@strata)
    if (!all(unlist(groups) %in% age_vec)) stop("Some group elements are not present in the strata.")

    rel_exp_mat <- rel_exp_matrix(phyex_set)
    
    # Create a properly formatted data frame for melting
    rel_exp_df <- data.frame(
        age = rep(as.integer(rownames(rel_exp_mat)), times = ncol(rel_exp_mat)),
        stage = rep(colnames(rel_exp_mat), each = nrow(rel_exp_mat)),
        expr = as.vector(rel_exp_mat)
    )
    rel_exp_df$stage <- factor(rel_exp_df$stage, levels = as.character(phyex_set@conditions))

    pal <- PS_colours(phyex_set@num_strata)

    # Single group
    if (length(groups) == 1) {
        df <- dplyr::filter(rel_exp_df, age %in% groups[[1]])
        if (nrow(df) == 0) stop("No data found for the specified group. Check your group specification.")
        
        p <- ggplot(df, aes(x = stage, y = expr, group = factor(age), color = factor(age))) +
            geom_line(linewidth = 1.2, ...) +
            labs(x = phyex_set@conditions_label, y = "Relative Expression Level", color = "Age") +
            theme_minimal() +
            scale_color_manual(values = pal[groups[[1]]], name = "Age") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        if (!is.null(modules) && !is.null(modules$mid)) {
            p <- p + annotate("rect",
                xmin = min(modules$mid) - 0.5, xmax = max(modules$mid) + 0.5,
                ymin = -Inf, ymax = Inf, alpha = alpha, fill = "#4d004b")
        }
        return(p)
    }

    # Two groups
    df1 <- dplyr::filter(rel_exp_df, age %in% groups[[1]])
    df2 <- dplyr::filter(rel_exp_df, age %in% groups[[2]])
    
    if (nrow(df1) == 0) stop("No data found for group 1. Check your group specification.")
    if (nrow(df2) == 0) stop("No data found for group 2. Check your group specification.")
    
    ylims <- if (adjust_range) range(rel_exp_df$expr, na.rm = TRUE) else NULL

    create_plot <- function(df, group_idx, group_name) {
        ggplot(df, aes(x = stage, y = expr, group = factor(age), color = factor(age))) +
            geom_line(linewidth = 1.2, ...) +
            labs(x = phyex_set@conditions_label, y = "Relative Expression Level",
                 title = paste(phyex_set@name, group_name), color = "Age") +
            theme_minimal() +
            scale_color_manual(values = pal[groups[[group_idx]]], name = "Age") +
            scale_y_continuous(limits = ylims) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }

    p1 <- create_plot(df1, 1, "(Group 1)")
    p2 <- create_plot(df2, 2, "(Group 2)")

    if (!is.null(modules) && !is.null(modules$mid)) {
        module_rect <- annotate("rect",
            xmin = min(modules$mid) - 0.5, xmax = max(modules$mid) + 0.5,
            ymin = -Inf, ymax = Inf, alpha = alpha, fill = "#4d004b")
        p1 <- p1 + module_rect
        p2 <- p2 + module_rect
    }

    return(list(group1 = p1, group2 = p2))
}

#' @title Plot Mean Relative Expression Levels as Barplot
#' @description Plots mean relative expression levels for age category groups using a PhyloExpressionSet S7 object, with statistical testing.
#' @param phyex_set A PhyloExpressionSet object.
#' @param groups A list of integer vectors specifying age categories (e.g., phylostrata) for each group (2+ groups).
#' @param p_adjust_method P-value adjustment for multiple testing.
#' @param ... Further arguments passed to ggplot2 geoms.
#' @return ggplot2 object.
#' @export
plot_relative_expression_bar <- function(
    phyex_set,
    groups,
    p_adjust_method = NULL,
    ...) {
    
    if (is.null(groups) || !is.list(groups) || length(groups) < 2) stop("groups must be a list of at least 2 groups.")
    
    age_vec <- as.integer((phyex_set@strata))
    if (!all(unlist(groups) %in% age_vec)) stop("Some group elements are not present in the strata.")
    
    # Get relative expression matrix
    rel_exp_mat <- rel_exp_matrix(phyex_set)
    
    # Compute mean and SE for each group
    n_groups <- length(groups)
    mean_mat <- matrix(NA_real_, nrow = n_groups, ncol = ncol(rel_exp_mat))
    se_mat <- matrix(NA_real_, nrow = n_groups, ncol = ncol(rel_exp_mat))
    
    for (i in seq_along(groups)) {
        idx <- which(rownames(rel_exp_mat) %in% as.character(groups[[i]]))
        vals <- rel_exp_mat[idx, , drop = FALSE]
        mean_mat[i, ] <- colMeans(vals, na.rm = TRUE)
        se_mat[i, ] <- apply(vals, 2, function(x) sd(x, na.rm = TRUE) / sqrt(length(x)))
    }
    
    # Kruskal-Wallis test per condition
    pvals <- vapply(seq_len(ncol(rel_exp_mat)), function(j) {
        group_vals <- lapply(groups, function(g) as.numeric(rel_exp_mat[as.character(g), j]))
        group_vals <- lapply(group_vals, function(x) x[is.finite(x)])
        
        if (all(lengths(group_vals) > 1)) {
            tryCatch(stats::kruskal.test(group_vals)$p.value, error = function(e) NA_real_)
        } else {
            NA_real_
        }
    }, numeric(1))
    
    if (!is.null(p_adjust_method)) pvals <- p.adjust(pvals, method = p_adjust_method)
    
    # Convert p-values to significance stars
    pval_stars <- case_when(
        is.na(pvals) ~ "",
        pvals > 0.05 ~ "NS",
        pvals <= 0.0005 ~ "***",
        pvals <= 0.005 ~ "**",
        pvals <= 0.05 ~ "*",
        .default = ""
    )
    
    # Prepare data for plotting
    df_bar <- data.frame(
        stage = rep(as.character(phyex_set@conditions), each = n_groups),
        mean = as.vector(mean_mat),
        se = as.vector(se_mat),
        group = factor(rep(seq_along(groups), times = phyex_set@num_conditions))
    )
    
    # Use PS_colours for consistent palette
    pal <- PS_colours(phyex_set@num_strata)
    
    # Create plot
    p <- ggplot(df_bar, aes(x = stage, y = mean, fill = group)) +
        geom_bar(stat = "identity", position = position_dodge(width = 0.7), 
                 color = "black", width = 0.7, ...) +
        geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                     width = 0.25, position = position_dodge(width = 0.7)) +
        labs(x = phyex_set@conditions_label, y = "Mean Relative Expression", 
             title = phyex_set@name, fill = "Age Group") +
        scale_fill_manual(values = pal[seq_len(n_groups)], 
                         labels = sapply(groups, function(g) paste0(min(g), "-", max(g)))) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "top"
        )
    
    # Add significance annotations
    if (any(pval_stars != "")) {
        p <- p + geom_text(
            data = data.frame(
                stage = as.character(phyex_set@conditions), 
                y = apply(mean_mat + se_mat, 2, max, na.rm = TRUE) + 0.02,
                label = pval_stars
            ),
            aes(x = stage, y = y, label = label), 
            inherit.aes = FALSE, vjust = 0, size = 3
        )
    }
    
    return(p)
}

