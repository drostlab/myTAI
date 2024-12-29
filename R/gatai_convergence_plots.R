
remove_genes <- function(phyex_set,
                         genes) {
    genes <- as.character(genes)
    filtered_set <- dplyr::filter(phyex_set, !GeneID %in% genes)
    return(filtered_set)
}

consensus <- function(x, p=0.5) {
    all <- Reduce(union, x)
    k <- length(x)
    counts <- all |>
        purrr::map_int(\(i) sum(purrr::map_lgl(x, \(v) i %in% v)))
    majority <- all[counts > k*p]
    
    return(majority)
}

full_gatai_convergence_plot <- function(phyex_set, 
                                        runs,
                                        p=0.5,
                                        ps=c(0.25,0.5,0.75)) {
    
    c <- convergence_plots(phyex_set, runs, ps)
    t <- threshold_comparison_plots(phyex_set, runs)
    petal <- petal_plot(runs, p=p)
    
    legend <- cowplot::get_legend(c[[1]])
    
    p <- cowplot::plot_grid(cowplot::plot_grid(c[[1]] + ggplot2::theme(legend.position = "none") +
                                                   ggplot2::labs(title="Convergence of GATAI consensus set sizes for different thresholds"), 
                                               c[[2]] + ggplot2::theme(legend.position = "none") + 
                                                   ggplot2::labs(title="Convergence of GATAI consensus set p values for different thresholds"),
                                               ncol=1),
                            legend,
                            cowplot::plot_grid(t[[1]] + ggplot2::labs(title="How the threshold of GATAI consensus influences set size"),
                                               t[[2]] + ggplot2::labs(title="How the threshold of GATAI consensus influences set p value"),
                                               ncol=1),
                            ncol=3,
                            rel_widths = c(2.2, 1, 2.4)
                            )
    p <- cowplot::plot_grid(p, 
                            cowplot::plot_grid(NULL,
                                               petal + ggplot2::labs(title="How many genes get lost per run."),
                                               NULL,
                                               rel_widths=c(0.3,0.4,0.3), 
                                               ncol=3), 
                            ncol=1, 
                            rel_heights = c(0.75,0.2)
                            )
    return(p)
}

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
    running_count <- Vectorize(running_count, vectorize.args = c("i", "p"))
    
    count_matrix <- outer(1:num_runs, ps, running_count)
    colnames(count_matrix) <- labels
    
    # Generate running p value matrix
    
    flt_null_dist <- flt_null_dist(phyex_set)
    
    running_pval <- function(i, p) {
        genes <- consensus(runs[1:i], p)
        filtered_phyex <- remove_genes(phyex_set, genes)
        pval <- flt_p_val(filtered_phyex, flt_null_dist)
        
        return(log10(pval))
    }
    running_pval <- Vectorize(running_pval, vectorize.args = c("i", "p"))
    
    pval_matrix <- outer(1:num_runs, ps, running_pval)
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
            ggplot2::scale_y_continuous(limits=ylims)
        return(p)
    }
    
    
    counts_plot <- running_plot(count_matrix, "Number of genes", ylims=c(0, max(count_matrix)))
    pval_plot <- running_plot(pval_matrix, "log10(p_flt)", ylims=c(min(pval_matrix[, ncol(pval_matrix)]), 0))
    
    
    return(list(counts=counts_plot, pval=pval_plot))
}

threshold_comparison_plots <- function(phyex_set, runs) {
    num_runs <- length(runs)
    
    # Gene counts per threshold
    counts <- 1:num_runs |>
        sapply(\(i) length(consensus(runs, p=i/num_runs-0.00001)))
    
    flt_null_dist <- flt_null_dist(phyex_set)
    pval <- function(i) {
        genes <- consensus(runs, p=i/num_runs-0.00001)
        filtered_phyex <- remove_genes(phyex_set, genes)
        pval <- flt_p_val(filtered_phyex, flt_null_dist)
        
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
            ggplot2::scale_y_continuous(limits=ylims)
        return(p)
    }
    
    
    counts_plot <- threshold_plot(counts, "Number of genes", ylims=c(0, max(counts)))
    pval_plot <- threshold_plot(pvals, "log10(p_flt)", ylims=c(min(pvals), 0))
    
    return(list(counts=counts_plot, pval=pval_plot))
}


# 
# #organs <- c("brain", "cerebellum", "heart", "kidney", "liver", "ovary", "testis")
# organs <- c("brain")
# 
# for (o in organs) {
#     runs <- list()
#     num_runs <- 20
#     for (i in 1:num_runs) {
#         runs[[i]] <- readr::read_delim(sprintf("data/%s/extracted_genes_%s.txt", o, i), delim="\n", col_names=c("GeneID"))$GeneID
#     }
#     phylo_set <- read.csv(sprintf("data/%s/phylo_set.tsv", o), sep="\t")
#     convergence_plot(phylo_set, runs, ps=c(0.25, 0.5, 0.75))
#     ggsave(filename=sprintf("plots/%s convergence plot.png", o))
#     threshold=0.5
#     petal_plot(runs, p=threshold) + labs(title=sprintf("%s GATAI petal plot, threshold %s", o, threshold))
#     ggsave(filename=sprintf("plots/%s petal plot.png", o))
# }


