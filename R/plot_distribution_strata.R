
#' @title Plot Distribution of Genes Across Phylostrata
#' @description Create a bar plot showing the distribution of genes across phylostrata,
#' with options for showing observed vs. expected ratios.
#' 
#' @param strata Named factor vector of phylostratum assignments (names are gene IDs)
#' @param selected_gene_ids Character vector of gene IDs to include in the plot
#' (default: all genes in strata)
#' @param as_log_obs_exp Logical indicating whether to show log2(observed/expected) ratios
#' instead of raw counts (default: FALSE)
#' 
#' @return A ggplot2 object showing the phylostratum distribution
#' 
#' @details
#' This function visualizes how genes are distributed across different phylostrata.
#' When as_log_obs_exp=FALSE, it shows raw gene counts per stratum. When TRUE, it
#' shows log2 ratios of observed vs. expected gene counts, useful for identifying
#' enrichment or depletion of specific strata in gene sets.
#' 
#' @examples
#' # Plot raw gene counts by strata
#' # p1 <- plot_distribution_strata(phyex_set@strata)
#' 
#' # Plot observed vs expected ratios for selected genes
#' # p2 <- plot_distribution_strata(phyex_set@strata, 
#' #                                selected_gene_ids = my_genes,
#' #                                as_log_obs_exp = TRUE)
#' 
#' @import ggplot2 dplyr tidyr
#' @export
plot_distribution_strata <- function(strata,
                                     selected_gene_ids = names(strata),
                                     as_log_obs_exp = FALSE) {

    if (!is.factor(strata)) {
        stop("'strata' must be a factor vector.", call. = FALSE)
    }

    df <- data.frame(Stratum=strata, GeneID=names(strata))
    df_selected <- df |> filter(GeneID %in% selected_gene_ids)
    if (!as_log_obs_exp) {
        ggplot(df_selected, aes(Stratum, fill=Stratum)) + 
            geom_bar(stat="count", position = "dodge", show.legend = TRUE) + 
            labs(fill="Stratum") +
            xlab("Stratum") + 
            ylab("Gene Count") +
            scale_x_discrete(labels=labels(df$Stratum), drop = FALSE) +
            scale_fill_manual(values = PS_colours(length(levels(strata))), labels=levels(df_selected$Stratum), drop=FALSE) +
            guides(fill=guide_legend(override.aes = list(size = 1.0), keyheight = 0.5, keywidth = 0.5)) +
            theme_minimal() +
            theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
            geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size=3)
    }
    else {
        counts_all <- df |> count(Stratum) |> mutate(p_all = n / sum(n))
        counts_sel <- df_selected |> count(Stratum) |> mutate(p_sel = n / sum(n))
        log_ratio <- full_join(counts_all, counts_sel, by= "Stratum") |>
            mutate(log_obs_exp = log2((p_sel + 1e-6)/ (p_all + 1e-6)))
        
        ggplot(log_ratio, aes(Stratum, 
                              y=log_obs_exp,
                              fill= log_obs_exp)) +
            geom_bar(stat = "identity", size=0.1) +
            labs(fill="Stratum") +
            xlab("Stratum") +
            ylab("Log(Obs/Exp)") +
            scale_x_discrete(labels=labels(df$Stratum), drop = FALSE) +
            scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
            guides(fill="none") + 
            theme_minimal()
    
        
    }
}

#' @title Calculate Phylostratum Enrichment
#' @description Calculate log2(observed/expected) enrichment ratios for phylostrata
#' in a selected gene set compared to the background distribution.
#' 
#' @param strata Named factor vector of phylostratum assignments (names are gene IDs)
#' @param selected_gene_ids Character vector of gene IDs to test for enrichment
#' 
#' @return A data frame with columns:
#' \describe{
#'   \item{Stratum}{Phylostratum factor levels}
#'   \item{log_obs_exp}{Log2 ratio of observed vs expected proportions}
#' }
#' 
#' @details
#' This function calculates enrichment or depletion of phylostrata in a gene set
#' by comparing the observed proportion of each stratum in the selected genes
#' to the expected proportion based on the background distribution in all genes.
#' 
#' Positive values indicate enrichment (more genes than expected), while negative 
#' values indicate depletion (fewer genes than expected).
#' 
#' @examples
#' # Calculate enrichment for a gene set
#' # enrichment <- strata_enrichment(phyex_set@strata, my_gene_set)
#' # print(enrichment)
#' 
#' @export
strata_enrichment <- function(strata, selected_gene_ids) {
    df <- data.frame(Stratum=strata, GeneID=names(strata))
    df_selected <- df |> filter(GeneID %in% selected_gene_ids)
    counts_all <- df |> count(Stratum, .drop=F) |> mutate(p_all = n / sum(n))
    counts_sel <- df_selected |> count(Stratum, .drop=F) |> mutate(p_sel = n / sum(n))
    log_ratio <- full_join(counts_all, counts_sel, by= "Stratum") |>
        mutate(log_obs_exp = log2((p_sel)/ (p_all))) |>
        select(Stratum, log_obs_exp)
    return(log_ratio)
}