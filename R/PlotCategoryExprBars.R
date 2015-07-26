#' @title Plot the Expression Levels of each Age or Divergence Category as Barplot or Violinplot
#' @description This function visualizes the expression level distribution of each phylostratum during each time point or experiment
#' as barplot or violin plot enabling users to quantify the age (PS) or divergence category (DS) specific contribution to the
#' corresponding transcriptome.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param legendName a character string specifying whether "PS" or "DS" are used to compute relative expression profiles.
#' @param test.stat a logical value indicating whether a Benjamini-Hochberg adjusted \code{\link{kruskal.test}} should be applied to determine
#' significant differences in age or divergence category specific expression.
#' @param type type of age or divergence category comparison. Specifications can be \code{type = "category-centered"} or \code{type = "stage-centered"}.
#' @param distr.type format of visualizing age or divergence category specific expression distributions. Either \code{distr.type = "barplot"} or
#' \code{distr.type = "violin"}. 
#' @param log.expr a logical value specifying whether or not expression levels should be log2-transformed before visualization.
#' @author Hajk-Georg Drost
#' @details This way of visualizing the gene expression distribution of each age (PS) or divergence category (DS) during
#' all developmental stages or experiments allows users to detect specific age or divergence categories contributing significant
#' levels of gene expression to the underlying biological process (transcriptome).
#' 
#' This quantification allows to conclude that genes originating in specific PS or DS contribute significantly more to the overall transcriptome
#' than other genes originating from different PS or DS categories. More specialized analyses such as \code{\link{PlotMeans}}, \code{\link{PlotRE}},
#' \code{\link{PlotBarRE}}, etc. will then allow to study the exact mean expression patterns of these age or divergence categories.
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
#' \item \code{distr.type = "barplot"} This specification allows users to visualize the expression distribution of all PS or DS as barplot.
#' \item \code{distr.type = "violin"} This specification allows users to visualize the expression distribution of all PS or DS as violin plot.
#' }
#' @examples 
#' 
#' data(PhyloExpressionSetExample)
#' data(DivergenceExpressionSetExample)
#'
#' # category-centered visualization of PS specific expression level distributions (log-scale)
#' PlotCategoryExprBars(ExpressionSet = PhyloExpressionSetExample,
#'                      legendName    = "PS",
#'                      type          = "category-centered",
#'                      distr.type    = "barplot",
#'                      log.expr      = TRUE)
#'                      
#' 
#' \dontrun{
#' # stage-centered visualization of PS specific expression level distributions (log-scale)
#' PlotCategoryExprBars(ExpressionSet = PhyloExpressionSetExample,
#'                      legendName    = "PS",
#'                      distr.type    = "barplot",
#'                      type          = "stage-centered",
#'                      log.expr      = TRUE)
#' 
#'                      
#'                                                                
#' # category-centered visualization of PS specific expression level distributions (log-scale)
#' # as violoin plot
#' PlotCategoryExprBars(ExpressionSet = PhyloExpressionSetExample,
#'                      legendName    = "PS",
#'                      distr.type    = "violin",
#'                      type          = "stage-centered",
#'                      log.expr      = TRUE)
#'
#'
#'
#'
#' # analogous for DivergenceExpressionSets
#' PlotCategoryExprBars(ExpressionSet = DivergenceExpressionSetExample,
#'                      legendName    = "DS",
#'                      type          = "category-centered",
#'                      distr.type    = "barplot",
#'                      log.expr      = TRUE)
#'
#'}
#' @seealso \code{\link{PlotMeans}}, \code{\link{PlotRE}}, \code{\link{PlotBarRE}}, \code{\link{age.apply}}
#' @export         

PlotCategoryExprBars <- function(ExpressionSet,
                             legendName,
                             test.stat  = TRUE,
                             type       = "category-centered",
                             distr.type = "barplot",
                             log.expr   = FALSE){
        
        is.ExpressionSet(ExpressionSet)
        
        if (!is.element(legendName, c("PS","DS")))
                stop ("Please specify 'legendName' as either 'PS' or 'DS'.")
        
        if (!is.element(type, c("category-centered","stage-centered")))
                stop ("Please specify 'type' as either 'category-centered' or 'stage-centered'.")
        
        if (!is.element(distr.type, c("barplot","violin")))
                stop ("Please specify 'distr.type' as either 'barplot' or 'violin'.")
        
        ncols <- ncol(ExpressionSet)
        nPS <- length(names(table(ExpressionSet[ , 1])))
        
        # global variable definition
        PS <- DS <- value <- Stage <- NULL
        
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
                                p_stage.cetered <- p.adjust(as.numeric(age.apply(ExpressionSet, function(x) format(kruskal.test(data.frame(x))$p.value,digits = 3))), method = "BH") 
                        }
                        
                        pValNames <- rep("n.s.",ncols-2)
                        pValNames[which(p_stage.cetered <= 0.05)] <- "*"
                        pValNames[which(p_stage.cetered <= 0.005)] <- "**"
                        pValNames[which(p_stage.cetered <= 0.0005)] <- "***"
                        pValNames[which(is.na(pValNames))] <- "n.s."
                } 
                
        }
        
        colnames(ExpressionSet)[1] <- legendName
        # reshape ExpressionSet from wide-format to long-format
        
        if (legendName == "PS")
                ReshapedExpressionSet <- reshape2::melt(ExpressionSet[ ,c(1,3:ncols)], id.vars = "PS")
        
        else if (legendName == "DS")
                ReshapedExpressionSet <- reshape2::melt(ExpressionSet[ ,c(1,3:ncols)], id.vars = "DS")
        
        ReshapedExpressionSet[ , 1] <- factor(ReshapedExpressionSet[ , 1], ordered = TRUE)
        colnames(ReshapedExpressionSet)[2] <- "Stage"
       
        if (distr.type == "barplot"){
                
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
                                                ggplot2::facet_grid(. ~ PS, labeller = ggplot2::label_both)  + ggplot2::labs(x = "\nPhylostratum", y = "Expression Level\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) +
                                                ggplot2::theme(axis.title.x = ggplot2::element_text(angle = 90,vjust = 0.5, hjust = 1)) + ggplot2::theme_minimal() 
                                        
                                        
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
                                                ggplot2::facet_grid(. ~ PS, labeller = ggplot2::label_both)  + ggplot2::labs(x = "\nPhylostratum", y = "log2(expression level)\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal()
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
                                                ggplot2::facet_grid(. ~ DS, labeller = ggplot2::label_both)  + ggplot2::labs(x = "\nDivergence Stratum", y = "Expression Level\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal()
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
                                                ggplot2::facet_grid(. ~ DS, labeller = ggplot2::label_both)  + ggplot2::labs(x = "\nDivergence Stratum", y = "log2(expression level)\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal()
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
                                                ggplot2::facet_grid(. ~ PS, labeller = ggplot2::label_both)  + ggplot2::labs(x = "\nPhylostratum", y = "Expression Level\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal()
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
                                                ggplot2::facet_grid(. ~ PS, labeller = ggplot2::label_both)  + ggplot2::labs(x = "\nPhylostratum", y = "log2(expression level)\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal()
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
                                                ggplot2::facet_grid(. ~ DS, labeller = ggplot2::label_both)  + ggplot2::labs(x = "\nDivergence Stratum", y = "Expression Level\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal()
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
                                                ggplot2::facet_grid(. ~ DS, labeller = ggplot2::label_both)  + ggplot2::labs(x = "\nDivergence Stratum", y = "log2(expression level)\n") + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(3), color = I("black")) + 
                                                ggplot2::geom_point(stat = "summary", fun.y = "mean", size = I(2.2), color = I("orange")) + 
                                                ggplot2::theme(legend.position = "bottom") + ggplot2::theme_minimal()
                                }
                        } 
                }
        } 
        
        
        stat.result <- t(as.data.frame(pValNames))
        
        if (legendName == "PS") 
                colnames(stat.result) <- paste0("PS",1:nPS)
        
        if (legendName == "DS") 
                colnames(stat.result) <- paste0("DS",1:nPS)
        
        rownames(stat.result) <- type
        
        print(stat.result)
        return (res)
}





