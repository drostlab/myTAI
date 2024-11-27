library(devtools)

load_all()

library(readr)
ExpressionSet <- read_delim("embryogenesis.csv", delim = "\t", 
                                                 escape_double = FALSE, 
                                                 col_types = cols(Phylostratum = col_integer()), 
                                                 trim_ws = TRUE)

all.genes <- ExpressionSet[[2]]
top20.genes <- TopExpressionGenes(ExpressionSet, p=0.8)
top10.genes <- TopExpressionGenes(ExpressionSet, p=0.9)
top20var.genes <- TopVarianceGenes(ExpressionSet, p=0.8)
top10var.genes <- TopVarianceGenes(ExpressionSet, p=0.9)
bottom80.genes <- subset(ExpressionSet, !(ExpressionSet[[2]] %in% top20.genes))[[2]]
gatai.genes <- scan("embryogenesis.txt", sep="\n", what='character')

all.set <- ExpressionSet
bottom90.set <-subset(ExpressionSet, !(ExpressionSet[[2]] %in% top10.genes))
bottom80.set <-subset(ExpressionSet, (ExpressionSet[[2]] %in% bottom80.genes))
bottom90var.set <-subset(ExpressionSet, !(ExpressionSet[[2]] %in% top10var.genes))
bottom80var.set <-subset(ExpressionSet, !(ExpressionSet[[2]] %in% top20var.genes))
gatai.set <- subset(ExpressionSet, !(ExpressionSet[[2]] %in% gatai.genes))
    
PlotSelectedAgeDistr(ExpressionSet, all.genes, legendName = "PS")
PlotSelectedAgeDistr(ExpressionSet, bottom80.genes, legendName = "PS")
PlotSelectedAgeDistr(ExpressionSet, top20.genes, legendName = "PS")
PlotSelectedAgeDistr(ExpressionSet, top20var.genes, legendName = "PS")

length(top10var.genes)
length(intersect(top10var.genes, top10.genes))


library(ggplot2)
library(patchwork)
PlotPermutedTAI <- function(ExpressionSets,
                            set.labels,
                            permutations=1000
) {
    dfs <- list()
    for (i in 1:length(ExpressionSets)){
        ExpressionSet <- ExpressionSets[[i]]
        set.label <- set.labels[[i]]
        bm <- bootMatrix(ExpressionSet, permutations=permutations)
        vars <- apply(bm, 1, stats::var)
        means <- apply(bm, 1, mean)
        tai <- TAI(ExpressionSet)
        tai_var <- stats::var(tai)
        tai_mean <- mean(tai)
        
        df <- tibble::as_tibble(bm) |>
            #tibble::rowid_to_column("Id") |>
            tibble::add_column(Group = "null dist", Label=set.label, Var=vars, Mean=means)

        df <- dplyr::bind_rows(tibble::tibble(Group = "true TAI", !!!tai, Label=set.label, Var=tai_var, Mean=tai_mean), df)
        
        #df_long <- tidyr::pivot_longer(df, cols=names(tai), names_to = "Stage", values_to = "Value")
        
        dfs[[i]] <- df
    }
    df <- dfs |>
        dplyr::bind_rows() |>
        tibble::rowid_to_column("Id")
    print.data.frame(head(df))
    
    
    df_long <- tidyr::pivot_longer(df, cols=names(tai), names_to = "Stage", values_to = "Value")
    print.data.frame(head(df_long))
    #df_long <- dplyr::arrange(df_long, dplyr::desc(df_longId))
    
    profile_plot <-ggplot(data=df_long,
                          aes(group = interaction(Id, Group), 
                              colour=factor(Label, levels=unique(Label)), 
                              linewidth=Group, 
                              alpha=Group, 
                              linetype=Group)) +
        geom_line(aes(x= factor(Stage, levels=unique(Stage)),
                      y = Value)) +
        labs(x = "Ontogeny", 
             y ="TAI", 
             colour="Expression Set", 
             title="Flat line test null distribution comparison for different gene quantiles",
             linewidth="Sample Type", alpha="Sample Type", linetype="Sample Type") +
        scale_linewidth_manual(values=c("true TAI" = 6.0, "null dist" = 0.2)) +
        scale_alpha_manual(values=c("true TAI" = 1.0, "null dist" = 0.2)) +
        scale_linetype_manual(values=c("true TAI" = "81", "null dist" = "solid")) +
        theme_minimal() +
        ggplot2::theme(
            title            = ggplot2::element_text(face = "bold", size=18),
            legend.title     = ggplot2::element_text(face = "bold", size=18),
            legend.text      = ggtext::element_markdown(size=14),
            axis.title       = ggplot2::element_text(face = "bold", size=20),
            axis.text.y      = ggplot2::element_text(face = "bold", size=18),
            axis.text.x      = ggplot2::element_text(face = "bold", size=18),
            panel.background = ggplot2::element_blank(),
            strip.text.x     = ggplot2::element_text(
                colour         = "black",
                face           = "bold",
                size           = 18
            )
        )
    
    variance_plot <-ggplot(data=df, 
                           aes(x=Mean,
                               fill=factor(Label, levels=unique(Label)))) +
        geom_vline(data=subset(df, df[[2]] == "true TAI"), 
                   aes(xintercept=Mean, 
                       colour=factor(Label, levels=unique(Label))),
                   linewidth=3,
                   linetype="81",
                   show.legend=TRUE) +
        geom_histogram(aes(y=..density..), bins=200, alpha=0.4, position="identity", colour="black") +
        labs(x = "Mean", 
             y ="Density", 
             title="Flat line test null distribution comparison for different gene quantiles",
             fill="Expression Set",
             colour="Expression Set") +
        theme_minimal() +
        ggplot2::theme(
            title            = ggplot2::element_text(face = "bold", size=18),
            legend.title     = ggplot2::element_text(face = "bold", size=18),
            legend.text      = ggtext::element_markdown(size=14),
            axis.title       = ggplot2::element_text(face = "bold", size=20),
            axis.text.y      = ggplot2::element_text(face = "bold", size=18),
            axis.text.x      = ggplot2::element_text(face = "bold", size=18),
            panel.background = ggplot2::element_blank(),
            strip.text.x     = ggplot2::element_text(
                colour         = "black",
                face           = "bold",
                size           = 18
            )
        )
        
    
    return(list(profile_plot = profile_plot, variance_plot = variance_plot))
}

PlotPermutedTAI(list(all.set, bottom90.set, bottom80.set), c("all", "bottom 90%", "bottom 80%"))

PlotGeneProfiles <- function(ExpressionSet) {
    vars <- .RowVars(ExpressionSet[, -2])
    means <- rowMeans(ExpressionSet[, -2])
    df <- tibble::as_tibble(ExpressionSet) |>
        tibble::add_column(Var=vars, Mean=means)
    print.data.frame(head(df))
    
    means_plot <- ggplot(df, aes(x=Mean)) +
        geom_histogram(bins=200, colour="black") +
        labs(x = "Mean", 
             y ="Count", 
             title="Distribution of mean gene expression for the genes") +
        theme_minimal() +
        ggplot2::theme(
            title            = ggplot2::element_text(face = "bold", size=18),
            legend.title     = ggplot2::element_text(face = "bold", size=18),
            legend.text      = ggtext::element_markdown(size=14),
            axis.title       = ggplot2::element_text(face = "bold", size=20),
            axis.text.y      = ggplot2::element_text(face = "bold", size=18),
            axis.text.x      = ggplot2::element_text(face = "bold", size=18),
            panel.background = ggplot2::element_blank(),
            strip.text.x     = ggplot2::element_text(
                colour         = "black",
                face           = "bold",
                size           = 18
            )
        )
    
    return(list(means_plot=means_plot))
        
}



PlotGeneProfiles(all.set)


