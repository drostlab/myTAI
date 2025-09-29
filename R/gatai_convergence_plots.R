#' @title Calculate Consensus Gene Set
#' @description Calculate consensus genes from multiple GATAI runs based on
#' frequency threshold.
#' 
#' @param x List of gene sets from different GATAI runs
#' @param p Frequency threshold (default: 0.5)
#' 
#' @return Character vector of consensus genes
#' 
#' @details
#' This function identifies genes that appear in at least p proportion of
#' the input gene sets, providing a consensus set of genes across multiple
#' GATAI runs.
#' 
#' @examples
#' # Calculate consensus from multiple runs
#' # consensus_genes <- consensus(gatai_runs, p = 0.5)
#' 
#' @keywords internal
consensus <- function(x, p=0.5) {
    all <- Reduce(union, x)
    k <- length(x)
    counts <- sapply(all, function(gene) {
        sum(sapply(x, function(gene_set) gene %in% gene_set))
    })
    majority <- all[counts >= k*p]
    
    return(majority)
}

#' @title Create Full GATAI Convergence Plot
#' @description Create a comprehensive plot showing GATAI convergence across multiple
#' runs and thresholds.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param runs List of GATAI run results
#' @param p Consensus threshold for petal plot (default: 0.5)
#' @param ps Vector of consensus thresholds for convergence plots (default: c(0.25, 0.5, 0.75))
#' 
#' @return A patchwork composition showing convergence analysis
#' 
#' @details
#' This function creates a comprehensive visualization of GATAI convergence
#' including consensus set sizes, p-values, threshold comparisons, and
#' gene removal patterns across multiple runs.
#' 
#' @examples
#' # Create full convergence plot
#' # conv_plot <- full_gatai_convergence_plot(phyex_set, gatai_runs)
#' 
#' @import patchwork
#' @keywords internal
full_gatai_convergence_plot <- function(phyex_set, 
                                        runs,
                                        p=0.5,
                                        ps=c(0.25,0.5,0.75)) {
    
    c <- convergence_plots(phyex_set, runs, ps)
    t <- threshold_comparison_plots(phyex_set, runs)
    petal <- petal_plot(runs, p=p)
    
    # Create plots with titles
    c1_titled <- c[[1]] + 
        ggplot2::labs(title="Convergence of GATAI consensus set sizes for different thresholds") +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 8))
    
    c2_titled <- c[[2]] + 
        ggplot2::labs(title="Convergence of GATAI consensus set p values for different thresholds") +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 8), legend.position="none")
    
    t1_titled <- t[[1]] + 
        ggplot2::labs(title="How the threshold of GATAI consensus influences set size") +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 8))
    
    t2_titled <- t[[2]] + 
        ggplot2::labs(title="How the threshold of GATAI consensus influences set p value") +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 8))
    
    petal_titled <- petal + 
    ggplot2::labs(title = "How many genes get lost per run.") +
    ggplot2::theme(
        plot.title = ggplot2::element_text(size = 8),
        aspect.ratio = 0.6  # Makes plot taller relative to width
    )
    
    # Create layout using patchwork
    # Top row with three columns: convergence plots, threshold plots
    top_left <- c1_titled / c2_titled
    top_right <- t1_titled / t2_titled
    
    # Combine top row with desired widths
    top_row <- top_left | top_right
    top_row <- top_row + patchwork::plot_layout(widths = c(2.2, 2.4))
    
    # Create empty spacer plots for whitespace and centering
    spacer_vertical <- ggplot2::ggplot() + ggplot2::theme_void()
    spacer_left <- ggplot2::ggplot() + ggplot2::theme_void()
    spacer_right <- ggplot2::ggplot() + ggplot2::theme_void()
    
    # Create centered petal plot row with spacers on sides
    petal_row <- spacer_left | petal_titled | spacer_right
    petal_row <- petal_row + patchwork::plot_layout(widths = c(1, 2, 1))  # Center with 2:1:2 ratio
    
    # Full layout: top row, vertical spacer, then centered petal plot
    p <- top_row / spacer_vertical / petal_row + 
        patchwork::plot_layout(heights = c(0.6, 0.1, 0.3))
    return(p)
}

#' @title Create Petal Plot for Gene Removal Analysis
#' @description Create a petal plot showing how many genes are removed per run
#' relative to the consensus set.
#' 
#' @param sets List of gene sets from GATAI runs
#' @param p Consensus threshold (default: 0.5)
#' 
#' @return A ggplot2 petal plot
#' 
#' @details
#' This function creates a petal plot visualization showing the relationship
#' between individual GATAI runs and the consensus gene set, highlighting
#' how many genes are unique to each run.
#' 
#' @examples
#' # Create petal plot
#' # petal <- petal_plot(gatai_runs, p = 0.5)
#' 
#' @keywords internal
petal_plot <- function(sets, p=0.5) {
    ref_set <- consensus(sets, p=p)
    num_sets <- length(sets)
    set_sizes <- sapply(sets, length)
    intersection_sizes <- sapply(sets, \(x) length(intersect(x, ref_set)))
    ref_set_size <- length(ref_set)
    
    df <- tibble::tibble(
        Set = rep(c("Consensus Set", paste0("Set ", 1:num_sets)), each=2),
        Type = rep(c("Original Set", "Intersection"), num_sets+1),
        Size = c(ref_set_size, 0, rbind(set_sizes, intersection_sizes))
    )
    
    p <- ggplot2::ggplot(df, 
                         ggplot2::aes(x=factor(Set, levels=unique(Set)), 
                                      y=Size, 
                                      fill=interaction(Type, Set == "Consensus Set"))
                         ) +
        ggplot2::geom_bar(stat="identity", 
                          position=ggplot2::position_identity()) +
        ggplot2::labs(x="Sets", y="Count") +
        ggplot2::scale_fill_manual(
            values = c(
                "Original Set.FALSE" = "#66c2a5",  
                "Intersection.FALSE" = "#fc8f92",  
                "Original Set.TRUE" = "#fc8d62"   
            ),
            labels=c("Intersection Sizes", "Set Sizes", "Consensus Set Size"),
            name="Counts"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            panel.background = ggplot2::element_rect(fill = "white"),
            plot.background = ggplot2::element_rect(fill = "white"),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )
    
    return(p)
}

#' @title Create Convergence Plots for GATAI Analysis
#' @description Create plots showing how consensus gene set sizes and p-values
#' converge across GATAI runs for different threshold values.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param runs List of GATAI run results
#' @param ps Vector of consensus thresholds to analyze (default: c(0.5))
#' 
#' @return A list with two ggplot objects:
#' \describe{
#'   \item{counts}{Plot showing convergence of consensus set sizes}
#'   \item{pval}{Plot showing convergence of p-values}
#' }
#' 
#' @details
#' This function analyzes how consensus gene sets and their statistical
#' significance change as more GATAI runs are included in the analysis.
#' It uses cached null distributions for efficient p-value calculation.
#' 
#' @examples
#' # Create convergence plots
#' # conv_plots <- convergence_plots(phyex_set, gatai_runs, ps = c(0.25, 0.5, 0.75))
#' 
#' @import ggplot2 tidyr patchwork
#' @keywords internal
convergence_plots <- function(phyex_set, runs, ps=c(0.5)) {
    num_runs <- length(runs)
    
    labels <- paste0("p=", ps)
    
    # Add union and intersection
    ps <- c(0.0, ps, 0.99999)
    labels <- c("p=0.0 (Union)", labels, "p=1.0 (Intersection)")
    
    # Generate running count matrix
    running_count <- function(i, p) {
        return(length(consensus(runs[1:i], p)))
    }
    running_count <- Vectorize(running_count, vectorize.args = c("i"))
    
    count_matrix <- sapply(ps, function(p) {
        sapply(1:num_runs, function(i) running_count(i, p))
    })
    colnames(count_matrix) <- labels
    
    # Generate running p value matrix using cached null distribution
    running_pval <- function(i, p) {
        genes <- consensus(runs[1:i], p)
        filtered_phyex <- remove_genes(phyex_set, genes)
        
        # Use flatline test which has cached null distribution
        test_result <- stat_flatline_test(filtered_phyex, plot_result = FALSE)
        pval <- test_result@p_value
        
        return(log10(pval))
    }
    running_pval <- Vectorize(running_pval, vectorize.args = c("i"))
    
    pval_matrix <- sapply(ps, function(p) {
        sapply(1:num_runs, function(i) running_pval(i, p))
    })
    colnames(pval_matrix) <- labels
    
    # Generate the 2 plots
    running_plot <- function(val_matrix, y_label, ylims) {
        df <- as.data.frame(val_matrix)
        df$Run <- 1:num_runs
        
        df <- df |>
            tidyr::pivot_longer(
                cols=!Run,
                names_to="Threshold",
                values_to=y_label
            )
        df$Threshold <- factor(df$Threshold, levels=unique(df$Threshold))
        
        p <- ggplot2::ggplot(df, 
                            ggplot2::aes(x=Run, 
                                         y=.data[[y_label]], 
                                         color=Threshold)
                            ) +
            ggplot2::geom_line(linewidth=1.3) +
            ggplot2::geom_hline(yintercept = 0) +
            ggplot2::scale_y_continuous(limits=ylims) +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
        return(p)
    }
    
    
    counts_plot <- running_plot(count_matrix, "Number of genes", ylims=c(0, max(count_matrix)))
    pval_plot <- running_plot(pval_matrix, "log10(p_flt)", ylims=c(min(pval_matrix[, ncol(pval_matrix)]), 0))
    
    
    return(list(counts=counts_plot, pval=pval_plot))
}

#' @title Create Threshold Comparison Plots
#' @description Create plots comparing how different consensus thresholds
#' affect gene set sizes and p-values in GATAI analysis.
#' 
#' @param phyex_set A PhyloExpressionSet object
#' @param runs List of GATAI run results
#' 
#' @return A list with two ggplot objects:
#' \describe{
#'   \item{counts}{Plot showing gene counts across thresholds}
#'   \item{pval}{Plot showing p-values across thresholds}
#' }
#' 
#' @details
#' This function analyzes how the choice of consensus threshold (minimum
#' number of runs a gene must appear in) affects the final gene set size
#' and statistical significance. Uses cached null distributions for
#' efficient p-value calculation.
#' 
#' @examples
#' # Create threshold comparison plots
#' # thresh_plots <- threshold_comparison_plots(phyex_set, gatai_runs)
#' 
#' @import ggplot2 patchwork
#' @keywords internal
threshold_comparison_plots <- function(phyex_set, runs) {
    num_runs <- length(runs)
    
    # Gene counts per threshold
    counts <- 1:num_runs |>
        sapply(\(i) length(consensus(runs, p=i/num_runs-0.00001)))
    
    # P-values using cached null distribution
    pval <- function(i) {
        genes <- consensus(runs, p=i/num_runs-0.00001)
        filtered_phyex <- remove_genes(phyex_set, genes)
        
        # Use flatline test which has cached null distribution
        test_result <- stat_flatline_test(filtered_phyex, plot_result = FALSE)
        pval <- test_result@p_value
        
        return(log10(pval))
    }
    pvals <- sapply(1:num_runs, pval)
    
    # Generate the 2 plots
    threshold_plot <- function(val_vec, y_label, ylims) {
        
        df <- data.frame(Threshold=1:num_runs, y=val_vec)
        
        p <- ggplot2::ggplot(df, 
                             ggplot2::aes(x=Threshold, 
                                          y=y)
                             ) +
            ggplot2::geom_line(linewidth=1.3) +
            ggplot2::labs(x="Threshold (gene in at least x runs)",y=y_label) +
            ggplot2::geom_hline(yintercept = 0) + 
            ggplot2::scale_y_continuous(limits=ylims) +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
        return(p)
    }
    
    
    counts_plot <- threshold_plot(counts, "Number of genes", ylims=c(0, max(counts)))
    pval_plot <- threshold_plot(pvals, "log10(p_flt)", ylims=c(min(pvals), 0))
    
    return(list(counts=counts_plot, pval=pval_plot))
}
