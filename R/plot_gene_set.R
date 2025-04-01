#' @title Plot the Expression Profiles of a Gene Set
#' @description
#' This function simply visualizes the gene expression profiles of
#' a defined subset of genes stored in the input \code{ExpressionSet}. The maximum nummber of genes to visualize
#' at once is 25.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param gene_set a character vector storing the gene ids for which gene expression profiles shall be visualized. 
#' @param get_subset a logical value indicating whether or not an \code{ExpressionSet} subset of the selected \code{gene_set} should be retuned. 
#' @param use_only_map a logical value indicating whether instead of a standard \code{ExpressionSet} only a \code{Phylostratigraphic Map} or \code{Divergene Map} is passed to the function.
#' @param colors colors for gene expression profiles. Default: \code{colors = NULL}, hence default colours are used.
#' @param plot_legend a logical value indicating whether gene ids should be printed as legend next to the plot.
#' @param y_ticks a numeric value specifying the number of ticks to be drawn on the y-axis.
#' @param digits_ylab a numeric value specifying the number of digits shown for the expression levels on the y-axis.
#' @param add_expression_values a logical value indicating whether expression values should be displayed on the plot.
#' @param n_genes_for_distance a numeric value specifying the number of top genes to be selected based on the highest distance sum.
#' @param ... additional parameters passed to \code{\link{matplot}}.
#' @author Hajk-Georg Drost and Filipa Martins Costa
#' @details
#' 
#' This function simply visualizes or subsets the gene expression levels of a set of genes
#' that are stored in the input \code{ExpressionSet}.
#' 
#' @seealso \code{\link{SelectGeneSet}}, \code{\link{PlotEnrichment}}, \code{\link{DiffGenes}}  
#' @examples
#' data(PhyloExpressionSetExample)
#' 
#' # the best parameter setting to visualize this plot:
#' # png("test_png.png",700,400)
#' PlotGeneSet(ExpressionSet = PhyloExpressionSetExample, 
#'             gene_set      = PhyloExpressionSetExample[1:5, 2], 
#'             lty           = 1, 
#'             lwd           = 4,
#'             xlab          = "Ontogeny",
#'             ylab          = "Expression Level")
#' 
#' # dev.off()
#' 
#' # In case you would like to work with the expression levels
#' # of selected genes you can specify the 'get_subset' argument:
#' 
#' plot_gene_set(ExpressionSet = PhyloExpressionSetExample, 
#'             gene_set      = PhyloExpressionSetExample[1:5, 2], 
#'             get_subset    = TRUE)
#' 
#' 
#' # get a gene subset using only a phylostratihraphic map
#' ExamplePSMap <- PhyloExpressionSetExample[ , 1:2]
#' 
#' plot_gene_set(ExpressionSet = ExamplePSMap, 
#'             gene_set      = PhyloExpressionSetExample[1:5, 2], 
#'             get_subset    = TRUE,
#'             use_only_map  = TRUE)
#'             
#' # In case you want the expression values to appear for the 2 genes with a most differentiated profile
#'plot_gene_set(ExpressionSet = PhyloExpressionSetExample, 
#'             gene_set      = PhyloExpressionSetExample[1:5, 2], 
#'             add_expression_values =T,
#'             n_genes_for_distance = 2
#'             )
#' @export

plot_gene_set_old <- function(ExpressionSet, 
                          gene_set, 
                          get_subset   = FALSE,
                          use_only_map = FALSE,
                          colors       = NULL,
                          plot_legend  = TRUE,
                          y_ticks      = 6,
                          digits_ylab  = 1,
                          line_width = 1,
                          point_size = 2,
                          add_expression_values = FALSE,
                          n_genes_for_distance = 1,... ){
  
  ExpressionSet <- as.data.frame(ExpressionSet)
  
  if (length(gene_set) > 25) {
    warning("For visualization reasons, the maximum number of genes allowed at once is 25. The gene set has been truncated.", call. = FALSE)
    gene_set <- gene_set[1:25]  # trim gene set to the first 25 genes
  }
  
  if (!use_only_map)
    is.ExpressionSet(ExpressionSet)
  
  GeneSubSet_indixes <- stats::na.omit(match(tolower(gene_set), tolower(ExpressionSet[ , 2])))
  
  if (length(GeneSubSet_indixes) == 0)
    stop ("None of your input gene ids could be found in the ExpressionSet.", call. = FALSE)
  
  if (length(GeneSubSet_indixes) != length(gene_set))
    warning ("Only ",length(GeneSubSet_indixes), " out of your ", length(gene_set), " gene ids could be found in the ExpressionSet.", call. = FALSE)
  
  GeneSubSet <- ExpressionSet[GeneSubSet_indixes , ]
  ncols <- ncol(GeneSubSet)
  
  if(!is.null(colors)){
    if (length(colors) != length(gene_set))
      stop ("The number of colors and the number of genes do not match.", call. = FALSE)
  }
  
  # http://www.compbiome.com/2010/12/r-using-rcolorbrewer-to-colour-your.html
  if (is.null(colors))
    colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8,"Dark2"))(length(GeneSubSet_indixes))                   

  GeneSubSet_long <- reshape2::melt(GeneSubSet, 
                                      id.vars = names(GeneSubSet)[1:2], 
                                      measure.vars = colnames(GeneSubSet)[3:ncol(GeneSubSet)],
                                      variable.name = "Stage",                    
                                      value.name = "Expression"                   
    )
  
  # Function to calculate Euclidean distance between each point and all other points
  calculate_euclidean_distance_across_genes <- function(df) {
    # calculates pairwise Euclidean distances between all points, not just within GeneID
    df <- df |>
      dplyr::mutate(id = dplyr::row_number()) |>
      dplyr::mutate(
        distance_sum = purrr::map_dbl(id, function(i) {
          sum(sqrt((Expression[i] - Expression)^2 + (as.numeric(Stage[i]) - as.numeric(Stage))^2))
        })
      ) |>
      dplyr::select(-id)
  }
  
  # apply the Euclidean distance calculation across all GeneIDs
  GeneSubSet_long_with_density <- calculate_euclidean_distance_across_genes(GeneSubSet_long)
  
  # selecting the optimal point (least crowded across GeneIDs)
  optimal_points <- dplyr::ungroup(dplyr::slice_max(dplyr::group_by(GeneSubSet_long_with_density, GeneID), order_by = distance_sum, n = 1))
  print(optimal_points)
  
  GeneSubSet_long$Stage <- factor(GeneSubSet_long$Stage, levels =  colnames(GeneSubSet)[3:ncol(GeneSubSet)])

  p <- ggplot2::ggplot(GeneSubSet_long, ggplot2::aes(x = Stage, y = Expression, color = GeneID)) + 
    ggplot2::geom_line(ggplot2::aes(group = GeneID), linewidth = line_width) +        
    ggplot2::geom_point(size = point_size) +                              
    ggrepel::geom_text_repel( 
      data = optimal_points, 
      ggplot2::aes(label = GeneID),                        
      size = 3,                                   
      max.overlaps = 50,                          
      box.padding = 0.5,                          
      point.padding = 0.5,                        
      force = 2,                                  
      nudge_y = 0.3,                              
      direction = "both",                         
      segment.size = 0.2,                         
      segment.color = "grey50",                   
      show.legend = FALSE                         
    ) +            
    ggplot2::scale_color_manual(values = colors) +               
    ggplot2::labs(x = "Ontogeny", y = "Expression Level", color = "Genes") + 
    ggplot2::theme_minimal()+
    ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = y_ticks)) +
    ggplot2::scale_y_continuous(labels = scales::number_format(accuracy = 10^(-digits_ylab))) 
  
  # add legands on the side
  if (plot_legend) {
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(title = "Gene ID"))
  } else {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  # add expression values 
  if (add_expression_values) {
    top_N_genes <- optimal_points |>
      dplyr::arrange(dplyr::desc(distance_sum)) |>
      dplyr::slice_head(n = n_genes_for_distance)
    
    selected_gene_data <- dplyr::filter(GeneSubSet_long, is.element(GeneID, top_N_genes$GeneID))
    
    p <- p + ggplot2::geom_text(
      data = selected_gene_data, 
      ggplot2::aes(label = round(Expression, 2)), 
      vjust = 1.5, 
      size = 3,  
      show.legend = FALSE
    )
  }
  
  print(p)
  
  if (get_subset)
    return(GeneSubSet)

}
