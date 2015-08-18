#' @title Plot the Expression Levels of each Age or Divergence Category as Boxplot, Violinplot, or Dotplot
#' @description This function visualizes the expression level distribution of each phylostratum during each time point or experiment
#' as boxplot, dot plot, or violin plot enabling users to quantify the age (PS) or divergence (DS) category specific contribution to the
#' corresponding transcriptome.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param legendName a character string specifying whether "PS" or "DS" are used to compute relative expression profiles.
#' @param test.stat a logical value indicating whether a Benjamini-Hochberg adjusted \code{\link{kruskal.test}} should be applied to determine
#' significant differences in age or divergence category specific expression.
#' @param type type of age or divergence category comparison. Specifications can be \code{type = "category-centered"} or \code{type = "stage-centered"}.
#' @param distr.type format of visualizing age or divergence category specific expression distributions. Either \code{distr.type = "boxplot"}, \code{distr.type = "dotplot"}, or
#' \code{distr.type = "violin"}. 
#' @param log.expr a logical value specifying whether or not expression levels should internally be log2-transformed before visualization.
#' @param gene.set a character vector storing the gene ids for which gene expression levels shall be visualized.
#' @author Hajk-Georg Drost
#' @details This way of visualizing the gene expression distribution of each age (PS) or divergence (DS) category during
#' all developmental stages or experiments allows users to detect specific age or divergence categories contributing significant
#' levels of gene expression to the underlying biological process (transcriptome).
#' 
#' This quantification allows users to conclude that genes originating in specific PS or DS contribute significantly more to the overall transcriptome
#' than other genes originating from different PS or DS categories. More specialized analyses such as \code{\link{PlotMeans}}, \code{\link{PlotRE}},
#' \code{\link{PlotBarRE}}, etc. will then allow to study the exact mean expression patterns of these age or divergence categories.
#' 
#' The statistical quantification of differences between expression levels of different age or divergence categories
#' is done by performing a \code{\link{kruskal.test}} with Benjamini & Hochberg p-value adjustment for multiple comparisons.
#' 
#' \itemize{
#' \item \code{type = "category-centered"} Here, the \code{\link{kruskal.test}} quantifies the differences of gene expression between all combinations of age or divergence categories for each stage or experiment separately. Here, a significant p-value quantifies that there is at least one pairwise comparison for which age or divergence categories significantly differ in their gene expression distribution. This type of analysis allows users to detect stages or experiments that show high diviation between age or divergence category contributions to the overall transcriptome or no significant deviations of age or divergence categories, suggesting equal age or divergence category contributions to the overall transcriptome.   
#' \item \code{type = "stage-centered"} Here, the \code{\link{kruskal.test}} quantifies the differences of gene expression between all stages or experiments for each age or divergence category separately. Hence, the test quantifies whether or not the gene expression distribution of a single age or divergence category significantly changes throughout development or experiments. This type of analysis allows users to detect specific age or divergence categories that significantly change their expression levels throughout development or experiments.
#' }
#' 
#' Argument Specifications:
#' 
#' Argument: type
#' 
#' \itemize{
#' \item \code{type = "category-centered"} This specification allows users to compare the differences between all age or
#'  divergence categories during the same stage or experiment.
#'  \item \code{type = "stage-centered"} This specification allows users to compare the differences between all age or
#'  divergence categories between stages or experiments.
#' }
#' 
#'   
#' Argument: distr.type
#' 
#' \itemize{
#' \item \code{distr.type = "boxplot"} This specification allows users to visualize the expression distribution of all PS or DS as boxplot.
#' \item \code{distr.type = "violin"} This specification allows users to visualize the expression distribution of all PS or DS as violin plot.
#' \item \code{distr.type = "dotplot"} This specification allows users to visualize the expression distribution of all PS or DS as dot plot.
#' }
#' 
#' 
#' Finally, users can specify a \code{gene.set} (a subset of genes included in the input \code{ExpressioSet})
#' for which expression levels should be visualized as boxplot, dotplot, or violinplot.
#' 
#' @return 
#' 
#' A boxplot, violin plot, or dot plot visualizing the gene expression levels of
#' different PS or DS categories.
#' 
#' Furthermore, the statistical test results returned from the \code{\link{kruskal.test}} are printed
#' to the console.
#' 
#' (1) '*'   = P-Value <= 0.05 
#'
#' (2) '**'  = P-Value <= 0.005  
#'
#' (3) '***' = P-Value <= 0.0005  
#' 
#' (4) 'n.s.' = not significant = P-Value > 0.05
#' 
#' @examples 
#' 
#' data(PhyloExpressionSetExample)
#' data(DivergenceExpressionSetExample)
#'
#'\dontrun{
#'
#' # category-centered visualization of PS specific expression level distributions (log-scale)
#' PlotCategoryExpr(ExpressionSet = PhyloExpressionSetExample,
#'                      legendName    = "PS",
#'                      test.stat     = TRUE,
#'                      type          = "category-centered",
#'                      distr.type    = "boxplot",
#'                      log.expr      = TRUE)
#'                      
#' 
#' # stage-centered visualization of PS specific expression level distributions (log-scale)
#' PlotCategoryExpr(ExpressionSet = PhyloExpressionSetExample,
#'                      legendName    = "PS",
#'                      test.stat     = TRUE,
#'                      distr.type    = "boxplot",
#'                      type          = "stage-centered",
#'                      log.expr      = TRUE)
#' 
#'                      
#'                                                                
#' # category-centered visualization of PS specific expression level distributions (log-scale)
#' # as violoin plot
#' PlotCategoryExpr(ExpressionSet = PhyloExpressionSetExample,
#'                      legendName    = "PS",
#'                      test.stat     = TRUE,
#'                      distr.type    = "violin",
#'                      type          = "stage-centered",
#'                      log.expr      = TRUE)
#'
#'
#'
#'
#' # analogous for DivergenceExpressionSets
#' PlotCategoryExpr(ExpressionSet = DivergenceExpressionSetExample,
#'                      legendName    = "DS",
#'                      test.stat     = TRUE,
#'                      type          = "category-centered",
#'                      distr.type    = "boxplot",
#'                      log.expr      = TRUE)
#'
#'
#' # visualize the expression levels of 500 example genes
#' set.seed(234)
#' example.gene.set <- PhyloExpressionSetExample[sample(1:25260,500) , 2]
#' 
#' PlotCategoryExpr(ExpressionSet = PhyloExpressionSetExample,
#'                  legendName    = "PS",
#'                  test.stat     = TRUE,
#'                  type          = "category-centered",
#'                  distr.type    = "boxplot",
#'                  log.expr      = TRUE,
#'                  gene.set      = example.gene.set)
#'                  
#'}
#' @seealso \code{\link{PlotMeans}}, \code{\link{PlotRE}}, \code{\link{PlotBarRE}}, \code{\link{age.apply}},
#' \code{\link{pTAI}}, \code{\link{pTDI}}, \code{\link{pStrata}}, \code{\link{pMatrix}}, \code{\link{TAI}}, \code{\link{TDI}}
#' @export         

PlotCategoryExpr <- function(ExpressionSet,
                             legendName,
                             test.stat  = TRUE,
                             type       = "category-centered",
                             distr.type = "boxplot",
                             log.expr   = FALSE,
                             gene.set   = NULL){
        
        is.ExpressionSet(ExpressionSet)
        
        if (!is.element(legendName, c("PS","DS")))
                stop ("Please specify 'legendName' as either 'PS' or 'DS'.")
        
        if (!is.element(type, c("category-centered","stage-centered")))
                stop ("Please specify 'type' as either 'category-centered' or 'stage-centered'.")
        
        if (!is.element(distr.type, c("boxplot","violin","dotplot")))
                stop ("Please specify 'distr.type' as either 'boxplot', 'dotplot', or 'violin'.")
        
        if (!is.logical(log.expr))
                stop ("'log.expr' can only be TRUE or FALSE.")
        
        ncols <- ncol(ExpressionSet)
        nPS <- length(names(table(ExpressionSet[ , 1])))
        
        # global variable definition
        PS <- DS <- value <- Stage <- NULL
        
        # determine gene subset
        if (!is.null(gene.set)){
                GeneSubSet.indixes <- stats::na.omit(match(tolower(gene.set), tolower(ExpressionSet[ , 2])))
                
                if (length(GeneSubSet.indixes) == 0)
                        stop ("None of your input gene ids could be found in the ExpressionSet.")
                
                if (length(GeneSubSet.indixes) != length(gene.set))
                        warning ("Only ",length(GeneSubSet.indixes), " out of your ", length(gene.set), " gene ids could be found in the ExpressionSet.")
                
                GeneSubSet <- ExpressionSet[GeneSubSet.indixes , ]
                ncols <- ncol(GeneSubSet)
                
                ExpressionSet <- GeneSubSet
        }
        
        if (!log.expr)
                max.value <- max(ExpressionSet[ , 3:ncols])
        
        if (log.expr)
                max.value <- max(tf(ExpressionSet,log2)[ , 3:ncols])
        
        if (test.stat){
                
                if (type == "stage-centered"){
                        # perform a Kruskal Test to detect stages of significant PS or DS variation using BH adjusted p-values
                        
                        if (log.expr){
                                p_stage.cetered <- stats::p.adjust(as.numeric(age.apply(tf(ExpressionSet,log2), function(x) format(stats::kruskal.test(data.frame(x))$p.value,digits = 3))), method = "BH")       
                        }
                        
                        else if (!log.expr){
                                p_stage.cetered <- stats::p.adjust(as.numeric(age.apply(ExpressionSet, function(x) format(stats::kruskal.test(data.frame(x))$p.value,digits = 3))), method = "BH") 
                        }
                        
                        pValNames <- rep("n.s.",length(names(table(ExpressionSet[ , 1]))))
                        pValNames[which(p_stage.cetered <= 0.05)] <- "*"
                        pValNames[which(p_stage.cetered <= 0.005)] <- "**"
                        pValNames[which(p_stage.cetered <= 0.0005)] <- "***"
                        pValNames[which(is.na(pValNames))] <- "n.s."
                } 
                
                if (type == "category-centered"){
                        # perform a Kruskal Test to detect stages of significant PS or DS variation using BH adjusted p-values
                        
                        if (log.expr){
                                p_category.cetered <-   stats::p.adjust(as.numeric(apply(tf(ExpressionSet,log2)[ , 3:ncols], 2 , function(x) stats::kruskal.test(x, g = ExpressionSet[ , 1])$p.value)), method = "BH") 
                        }
                        
                        else if (!log.expr){
                                p_category.cetered <-  stats::p.adjust(as.numeric(apply(ExpressionSet[ , 3:ncols], 2 , function(x) stats::kruskal.test(x, g = ExpressionSet[ , 1])$p.value)), method = "BH") 
                        }
                        
                        pValNames <- rep("n.s.",ncols-2)
                        pValNames[which(p_category.cetered <= 0.05)] <- "*"
                        pValNames[which(p_category.cetered <= 0.005)] <- "**"
                        pValNames[which(p_category.cetered <= 0.0005)] <- "***"
                        pValNames[which(is.na(pValNames))] <- "n.s."
                }
        }
        
        colnames(ExpressionSet)[1] <- legendName
        # reshape ExpressionSet from wide-format to long-format
        
        if (legendName == "PS")
                ReshapedExpressionSet <- reshape2::melt(ExpressionSet[ ,c(1,3:ncols)], id.vars = "PS")
        
        else if (legendName == "DS")
                ReshapedExpressionSet <- reshape2::melt(ExpressionSet[ ,c(1,3:ncols)], id.vars = "DS")
        
        ReshapedExpressionSet[ , 1] <- factor(ReshapedExpressionSet[ , 1], ordered = TRUE, levels = 1:nPS)
        colnames(ReshapedExpressionSet)[2] <- "Stage"
       
        if (distr.type == "boxplot"){
                
                if (legendName == "PS"){
                        
                        if (!log.expr){
                                
                                if (type == "category-centered"){
                                        
                                        res <- ggplot2::ggplot(ReshapedExpressionSet, ggplot2::aes(x = PS, y = value, fill = PS))  + ggplot2::geom_boxplot() + 
                                                ggplot2::facet_grid(. ~ Stage, labeller = ggplot2::label_parsed)  + ggplot2::labs(x = "\nPhylostratum", y = "Expression Level\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal()
                                }
                                
                                else if (type == "stage-centered"){
                                        
#                                         pval_mapping <- data.frame(PS = 1:nPS, x_coord = rep(ifelse(nPS < 3, 1, 3),nPS),max_value = rep(max.value,nPS),pvals = pValNames)
#                                         print(pval_mapping)
                                        res <- ggplot2::ggplot(ReshapedExpressionSet, ggplot2::aes(x = Stage, y = value, fill = Stage))  + ggplot2::geom_boxplot() + 
                                                ggplot2::facet_grid(. ~ PS, labeller = ggplot2::label_both)  + ggplot2::labs(x = "\nDevelopmental Stage", y = "Expression Level\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) +
                                                ggplot2::theme(legend.position = "bottom") +
                                                ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1)) 
                                        
                                        
                                        # ggplot2::geom_text(data = pval_mapping, ggplot2::aes(x = x_coord, y = max_value,label = pvals), colour = "red", size = 5)
                                        
                                        # ggplot2::annotate("text",x = pval_mapping[ , "x_coord"], y = pval_mapping[ , "max_value"],label = pval_mapping[ , "pvals"], colour = "red", size = 5, fill = pval_mapping[ , "PS"])
                                        
                                        # ggplot2::geom_text(data = ReshapedExpressionSet, aes(label = pValNames, x = Stage, y = median(length(Stage))),  size=5)
                                                # ggplot2::geom_text(data = pval_mapping, ggplot2::aes(x = pval_mapping[ , "x_coord"], y = pval_mapping[ , "max_value"],label = pval_mapping[ , "pvals"]), colour = "red", size = 5)
                                                
                                        # ggplot2::annotate("text",x = 1:nPS, y = rep(max.value,nPS),label = pValNames, colour = "red", size = 5)
                                }
                        } 
                        
                        else if (log.expr){
                                
                                if (type == "category-centered"){
                                        
                                        res <- ggplot2::ggplot(ReshapedExpressionSet, ggplot2::aes(x = PS, y = log2(value), fill = PS))  + ggplot2::geom_boxplot() + 
                                                ggplot2::facet_grid(. ~ Stage, labeller = ggplot2::label_parsed)  + ggplot2::labs(x = "\nPhylostratum", y = "log2(expression level)\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal()
                                }
                                
                                else if (type == "stage-centered"){
                                        
                                        res <- ggplot2::ggplot(ReshapedExpressionSet, ggplot2::aes(x = Stage, y = log2(value), fill = Stage))  + ggplot2::geom_boxplot() + 
                                                ggplot2::facet_grid(. ~ PS, labeller = ggplot2::label_both)  + ggplot2::labs(x = "\nDevelopmental Stage", y = "log2(expression level)\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1)) 
                                }
                        } 
                }
                
                if (legendName == "DS"){
                        
                        if (!log.expr){
                                
                                if (type == "category-centered"){
                                        
                                        res <- ggplot2::ggplot(ReshapedExpressionSet, ggplot2::aes(x = DS, y = value, fill = DS))  + ggplot2::geom_boxplot() + 
                                                ggplot2::facet_grid(. ~ Stage, labeller = ggplot2::label_parsed)  + ggplot2::labs(x = "\nDivergence Stratum", y = "Expression Level\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal()
                                }
                                
                                else if (type == "stage-centered"){
                                        
                                        res <- ggplot2::ggplot(ReshapedExpressionSet, ggplot2::aes(x = Stage, y = value, fill = Stage))  + ggplot2::geom_boxplot() + 
                                                ggplot2::facet_grid(. ~ DS, labeller = ggplot2::label_both)  + ggplot2::labs(x = "\nDevelopmental Stage", y = "Expression Level\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1))
                                }
                        } 
                        
                        else if (log.expr){
                                
                                if (type == "category-centered"){
                                        
                                        res <- ggplot2::ggplot(ReshapedExpressionSet, ggplot2::aes(x = DS, y = log2(value), fill = DS))  + ggplot2::geom_boxplot() + 
                                                ggplot2::facet_grid(. ~ Stage, labeller = ggplot2::label_parsed)  + ggplot2::labs(x = "\nDivergence Stratum", y = "log2(expression level)\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal()
                                }
                                
                                else if (type == "stage-centered"){
                                        
                                        res <- ggplot2::ggplot(ReshapedExpressionSet, ggplot2::aes(x = Stage, y = log2(value), fill = Stage))  + ggplot2::geom_boxplot() + 
                                                ggplot2::facet_grid(. ~ DS, labeller = ggplot2::label_both)  + ggplot2::labs(x = "\nDevelopmental Stage", y = "log2(expression level)\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1))
                                }
                        } 
                }
        } 
        
        if (distr.type == "violin"){
                
                if (legendName == "PS"){
                        
                        if (!log.expr){
                                
                                if (type == "category-centered"){
                                        
                                        res <- ggplot2::ggplot(ReshapedExpressionSet, ggplot2::aes(x = PS, y = value, fill = PS))  + ggplot2::geom_violin(trim = FALSE) + 
                                                ggplot2::facet_grid(. ~ Stage, labeller = ggplot2::label_parsed)  + ggplot2::labs(x = "\nPhylostratum", y = "Expression Level\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal()
                                }
                                
                                else if (type == "stage-centered"){
                                        
                                        res <- ggplot2::ggplot(ReshapedExpressionSet, ggplot2::aes(x = Stage, y = value, fill = Stage))  + ggplot2::geom_violin(trim = FALSE) + 
                                                ggplot2::facet_grid(. ~ PS, labeller = ggplot2::label_both)  + ggplot2::labs(x = "\nDevelopmental Stage", y = "Expression Level\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1))
                                }
                        } 
                        
                        else if (log.expr){
                                
                                if (type == "category-centered"){
                                        
                                        res <- ggplot2::ggplot(ReshapedExpressionSet, ggplot2::aes(x = PS, y = log2(value), fill = PS))  + ggplot2::geom_violin(trim = FALSE) + 
                                                ggplot2::facet_grid(. ~ Stage, labeller = ggplot2::label_parsed)  + ggplot2::labs(x = "\nPhylostratum", y = "log2(expression level)\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal()
                                }
                                
                                else if (type == "stage-centered"){
                                        
                                        res <- ggplot2::ggplot(ReshapedExpressionSet, ggplot2::aes(x = Stage, y = log2(value), fill = Stage))  + ggplot2::geom_violin(trim = FALSE) + 
                                                ggplot2::facet_grid(. ~ PS, labeller = ggplot2::label_both)  + ggplot2::labs(x = "\nDevelopmental Stage", y = "log2(expression level)\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1))
                                }
                        } 
                }
                
                if (legendName == "DS"){
                        
                        if (!log.expr){
                                
                                if (type == "category-centered"){
                                        
                                        res <- ggplot2::ggplot(ReshapedExpressionSet, ggplot2::aes(x = DS, y = value, fill = DS))  + ggplot2::geom_violin(trim = FALSE) + 
                                                ggplot2::facet_grid(. ~ Stage, labeller = ggplot2::label_parsed)  + ggplot2::labs(x = "\nDivergence Stratum", y = "Expression Level\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal()
                                }
                                
                                else if (type == "stage-centered"){
                                        
                                        res <- ggplot2::ggplot(ReshapedExpressionSet, ggplot2::aes(x = Stage, y = value, fill = Stage))  + ggplot2::geom_violin(trim = FALSE) + 
                                                ggplot2::facet_grid(. ~ DS, labeller = ggplot2::label_both)  + ggplot2::labs(x = "\nDevelopmental Stage", y = "Expression Level\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1))
                                }
                        } 
                        
                        else if (log.expr){
                                
                                if (type == "category-centered"){
                                        
                                        res <- ggplot2::ggplot(ReshapedExpressionSet, ggplot2::aes(x = DS, y = log2(value), fill = DS))  + ggplot2::geom_violin(trim = FALSE) + 
                                                ggplot2::facet_grid(. ~ Stage, labeller = ggplot2::label_parsed)  + ggplot2::labs(x = "\nDivergence Stratum", y = "log2(expression level)\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal()
                                }
                                
                                else if (type == "stage-centered"){
                                        
                                        res <- ggplot2::ggplot(ReshapedExpressionSet, ggplot2::aes(x = Stage, y = log2(value), fill = Stage))  + ggplot2::geom_violin(trim = FALSE) + 
                                                ggplot2::facet_grid(. ~ DS, labeller = ggplot2::label_both)  + ggplot2::labs(x = "\nDevelopmental Stage", y = "log2(expression level)\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1))
                                }
                        } 
                }
        } 
        
        if (distr.type == "dotplot"){
                
                if (legendName == "PS"){
                        
                        if (!log.expr){
                                
                                if (type == "category-centered"){
                                        
                                        res <- ggplot2::ggplot(ReshapedExpressionSet, ggplot2::aes(x = PS, y = value, fill = PS))  + ggplot2::geom_dotplot(trim = FALSE) + 
                                                ggplot2::facet_grid(. ~ Stage, labeller = ggplot2::label_parsed)  + ggplot2::labs(x = "\nPhylostratum", y = "Expression Level\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal()
                                }
                                
                                else if (type == "stage-centered"){
                                        
                                        res <- ggplot2::ggplot(ReshapedExpressionSet, ggplot2::aes(x = Stage, y = value, fill = Stage))  + ggplot2::geom_dotplot(trim = FALSE) + 
                                                ggplot2::facet_grid(. ~ PS, labeller = ggplot2::label_both)  + ggplot2::labs(x = "\nDevelopmental Stage", y = "Expression Level\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1))
                                }
                        } 
                        
                        else if (log.expr){
                                
                                if (type == "category-centered"){
                                        
                                        res <- ggplot2::ggplot(ReshapedExpressionSet, ggplot2::aes(x = PS, y = log2(value), fill = PS))  + ggplot2::geom_dotplot(trim = FALSE) + 
                                                ggplot2::facet_grid(. ~ Stage, labeller = ggplot2::label_parsed)  + ggplot2::labs(x = "\nPhylostratum", y = "log2(expression level)\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal()
                                }
                                
                                else if (type == "stage-centered"){
                                        
                                        res <- ggplot2::ggplot(ReshapedExpressionSet, ggplot2::aes(x = Stage, y = log2(value), fill = Stage))  + ggplot2::geom_dotplot(trim = FALSE) + 
                                                ggplot2::facet_grid(. ~ PS, labeller = ggplot2::label_both)  + ggplot2::labs(x = "\nDevelopmental Stage", y = "log2(expression level)\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1))
                                }
                        } 
                }
                
                if (legendName == "DS"){
                        
                        if (!log.expr){
                                
                                if (type == "category-centered"){
                                        
                                        res <- ggplot2::ggplot(ReshapedExpressionSet, ggplot2::aes(x = DS, y = value, fill = DS))  + ggplot2::geom_dotplot(trim = FALSE) + 
                                                ggplot2::facet_grid(. ~ Stage, labeller = ggplot2::label_parsed)  + ggplot2::labs(x = "\nDivergence Stratum", y = "Expression Level\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal()
                                }
                                
                                else if (type == "stage-centered"){
                                        
                                        res <- ggplot2::ggplot(ReshapedExpressionSet, ggplot2::aes(x = Stage, y = value, fill = Stage))  + ggplot2::geom_dotplot(trim = FALSE) + 
                                                ggplot2::facet_grid(. ~ DS, labeller = ggplot2::label_both)  + ggplot2::labs(x = "\nDevelopmental Stage", y = "Expression Level\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1))
                                }
                        } 
                        
                        else if (log.expr){
                                
                                if (type == "category-centered"){
                                        
                                        res <- ggplot2::ggplot(ReshapedExpressionSet, ggplot2::aes(x = DS, y = log2(value), fill = DS))  + ggplot2::geom_dotplot(trim = FALSE) + 
                                                ggplot2::facet_grid(. ~ Stage, labeller = ggplot2::label_parsed)  + ggplot2::labs(x = "\nDivergence Stratum", y = "log2(expression level)\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal()
                                }
                                
                                else if (type == "stage-centered"){
                                        
                                        res <- ggplot2::ggplot(ReshapedExpressionSet, ggplot2::aes(x = Stage, y = log2(value), fill = Stage))  + ggplot2::geom_dotplot(trim = FALSE) + 
                                                ggplot2::facet_grid(. ~ DS, labeller = ggplot2::label_both)  + ggplot2::labs(x = "\nDevelopmental Stage", y = "log2(expression level)\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1))
                                }
                        } 
                }
        } 
        
        if (test.stat){
                stat.result <- t(as.data.frame(pValNames))
                
                if (legendName == "PS"){
                        if (type == "stage-centered"){
                                
                                if (is.null(gene.set))
                                        colnames(stat.result) <- paste0("PS",1:nPS)
                                
                                if (!is.null(gene.set))
                                        colnames(stat.result) <- paste0("PS",names(table(ExpressionSet[ , 1])))
                        }
                        
                        if (type == "category-centered")
                                colnames(stat.result) <- names(ExpressionSet)[3:ncols]
                } 
                
                if (legendName == "DS"){
                        if (type == "stage-centered"){
                                if (is.null(gene.set))
                                        colnames(stat.result) <- paste0("DS",1:nPS)
                                
                                if (!is.null(gene.set))
                                        colnames(stat.result) <- paste0("DS",names(table(ExpressionSet[ , 1])))
                        }
                        
                        if (type == "category-centered")
                                colnames(stat.result) <- names(ExpressionSet)[3:ncols]
                } 
                
                rownames(stat.result) <- type
                
                print(stat.result)
                
                if (!is.null(gene.set)){
                        cat("\n")
                        cat("# Genes: ")
                        print(table(GeneSubSet[ , 1]))
                }
        }
        
        return (res)
}





