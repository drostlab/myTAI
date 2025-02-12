

generate_null_txi_sample <- function(phyex_set,
                                     sample_size=10000) {
    age_vector <- as.vector(phyex_set[[1]])
    count_matrix <- as.matrix(phyex_set[3:ncol(phyex_set)])
    
    # columns: stages, rows: permuted instances, entries: TAI value
    permuted_txi_matrix <- cpp_bootMatrix(count_matrix,
                                          age_vector,
                                          sample_size)
    
    colnames(permuted_txi_matrix) <- colnames(count_matrix)
    
    return(permuted_txi_matrix)
}

plot_null_txi_sample <- function(null_txis, 
                                 test_txis) {
    S <- ncol(null_txis)
    if (unique(sapply(test_txis, length)) != S)
        stop("The number of stages between the null_txis and the test_txis is inconsistent")
    
    null_df <- tibble::as_tibble(null_txis) |>
        tibble::rowid_to_column("Id") |>
        tidyr::pivot_longer(cols=-Id, names_to = "Stage", values_to = "TXI") |>
        tibble::add_column(Group = "Null Hypothesis Sample")
    
    test_df <- tibble::tibble(
        Id = rep(1:length(test_txis), each = S),
        Group = rep(names(test_txis), each = S),
        Stage = rep(colnames(null_txis), times=length(test_txis)),
        TXI = unlist(test_txis)
    )
    
    avg <- mean(null_txis)
    
    colour_values <- setNames(RColorBrewer::brewer.pal(length(unique(test_df$Group)), "Set2"), unique(test_df$Group))
    colour_values["Null Hypothesis Sample"] <- "paleturquoise3"
    
    p <- ggplot2::ggplot(test_df,
                         ggplot2::aes(x=factor(Stage, levels=unique(Stage)), 
                                      y=TXI, 
                                      colour=factor(Group, levels=unique(Group)), 
                                      group=interaction(Group, Id))) +
        ggplot2::geom_line(linewidth=1.6) + 
        ggplot2::geom_line(data=null_df,
                           ggplot2::aes(x=factor(Stage, levels=unique(Stage)), 
                                        y=TXI,  
                                        group=interaction(Group, Id),
                                        colour=factor(Group, levels=unique(Group))),
                           linewidth=0.1, alpha=0.1) + 
        ggplot2::geom_hline(yintercept=avg, colour="paleturquoise4", linewidth=0.7) +
        ggplot2::scale_colour_manual(values=colour_values) +
        ggplot2::labs(x="Ontogeny",
                      y="TAI",
                      colour="Sample Type") + 
        ggplot2::scale_x_discrete(expand=c(0, 0)) +
        ggplot2::theme_minimal()
    
    return(p)
    
}