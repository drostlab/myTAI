#' @title Plot Relative Expression Levels
#' @description 
#' This function computes for each age category the corresponding relative expression profile.
#' 
#' For each age category the corresponding relative expression profile is being computed as follows:
#' 
#' \deqn{f_js = ( e_js - e_j min ) / ( e_j max - e_j min )}
#'
#' where \eqn{e_j min} and \eqn{e_j max} denote the minimum/maximum \code{\link{mean}} expression level 
#' of phylostratum j over  developmental stages s. This linear transformation corresponds to 
#' a shift by \eqn{e_j min} and a subsequent shrinkage by \eqn{e_j max - e_j min}. 
#' As a result, the relative expression level \eqn{f_js} of developmental stage s 
#' with minimum \eqn{e_js} is 0, the relative expression level \eqn{f_js} of the developmental 
#' stage s with maximum \eqn{e_js} is 1, and the relative expression levels \eqn{f_js} of 
#' all other stages s range between 0 and 1, accordingly.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param Groups a list containing the age categories for which mean expression levels shall be drawn.
#' For ex. evolutionary users can compare old phylostrata: PS1-3 (Class 1) and evolutionary young phylostrata: PS4-12 (Class 2). 
#' In this example, the list could be assigned as, \code{Groups = list(c(1:3), c(4:12))}. 
#' The group options is limited to 2 Groups.
#' @param modules a list storing three elements for specifying the modules: early, mid, and late. 
#' Each element expects a numeric vector specifying the developmental stages 
#' or experiments that correspond to each module. For example, 
#' \code{module} = \code{list(early = 1:2, mid = 3:5, late = 6:7)} devides a dataset storing seven developmental stages into 3 modules. Default is \code{modules = NULL}. 
#' But if specified, a shaded are will be drawn to illustrate stages corresponding to the mid module.
#' @param legendName a character string specifying the legend title.
#' @param xlab label of x-axis.
#' @param ylab label of y-axis.
#' @param main main text.
#' @param y.ticks number of ticks that shall be drawn on the y-axis.
#' @param adjust.range logical indicating whether or not the y-axis scale shall be adjusted to the same range in case two groups are specified. Default is \code{adjust.range = TRUE}.
#' @param alpha transparency of the shaded area (between [0,1]). Default is \code{alpha = 0.1}.
#' @param ... place holder for old version of PlotRE that was based on base graphics instead of ggplot2.
#' @details Studying the relative expression profiles of each phylostratum or divergence-stratum enables the detection
#' of common gene expression patterns shared by several phylostrata or divergence-strata.
#'
#' Finding similar relative expression profiles among phylostrata or divergence-strata suggests that 
#' phylostrata or divergence-strata sharing a similar relative expression profile are regulated by similar
#' gene regulatory elements. Hence, these common phylostrata or divergence-strata might govern similar processes in the given developmental time course. 
#' @return a plot showing the relative expression profiles of phylostrata or divergence-strata belonging to the same group.
#' @references 
#' Domazet-Loso T and Tautz D. 2010. "A phylogenetically based transcriptome age index mirrors ontogenetic divergence patterns". Nature (468): 815-818.
#'
#' Quint M et al. 2012. "A transcriptomic hourglass in plant embryogenesis". Nature (490): 98-101.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{PlotBarRE}}, \code{\link{RE}}, \code{\link{REMatrix}}
#' @examples
#' 
#' # read standard phylotranscriptomics data
#' data(PhyloExpressionSetExample)
#' data(DivergenceExpressionSetExample)
#'
#' # example PhyloExpressionSet
#' PlotRE(PhyloExpressionSetExample,
#'        Groups = list(c(1:3), c(4:12)), 
#'        legendName = "PS")
#'
#'
#' # or you can choose any combination of groups
#' PlotRE(PhyloExpressionSetExample,
#'        Groups = list(c(1,7,9), c(2:6,8,10:12)),
#'        legendName = "PS")
#' 
#'     
#' # example DivergenceExpressionSet
#' PlotRE(DivergenceExpressionSetExample,
#'        Groups = list(c(1:5), c(6:10)), 
#'        legendName = "DS")
#'
#'
#'   
#' @export

PlotRE <- function(ExpressionSet,
                   Groups     = NULL,
                   modules    = NULL,
                   legendName = "age",
                   xlab = "Ontogeny",
                   ylab = "Relative Expression Level",
                   main = "",
                   y.ticks = 10,
                   adjust.range = TRUE,
                   alpha = 0.008, ...)
{
        
        ExpressionSet <- as.data.frame(ExpressionSet)
        is.ExpressionSet(ExpressionSet)
        
        stage <- expr <- age <- NULL
        
        if(is.null(Groups))
                stop("Your Groups list does not store any items.", call. = FALSE)
        
        ### getting the PS names available in the given expression set
        age_names <- as.character(names(table(ExpressionSet[ , 1])))
        
        # test whether all group elements are available in the age vector
        # ra <- range(ExpressionSet[ , 1])
        if(!all(unlist(Groups) %in% as.numeric(age_names)))
                stop("There are items in your Group elements that are not available in the age column of your ExpressionSet.", call. = FALSE)
        
        if (length(Groups) > 2)
                stop("Please specify at maximum 2 groups that shall be compared.", call. = FALSE)
        
        ### getting the PS names available in the given expression set
        nPS <- length(age_names)
        nCols <- dim(ExpressionSet)[2]
        ### define and label the REmatrix that holds the rel. exp. profiles
        ### for the available PS
        MeanValsMatrix <- matrix(NA_real_,nPS,nCols-2)
        rownames(MeanValsMatrix) <- age_names
        colnames(MeanValsMatrix) <- names(ExpressionSet)[3:nCols]
        
        MeanValsMatrix <- age.apply(ExpressionSet, RE)
        mean.age <- data.frame(age = age_names, MeanValsMatrix, stringsAsFactors = FALSE)
        mMatrix <- tibble::as_tibble(reshape2::melt(mean.age, id.vars = "age"))
        colnames(mMatrix)[2:3] <- c("stage", "expr")
        
        if (length(Groups) == 1) {
                
                p <- ggplot2::ggplot(mMatrix, ggplot2::aes( factor(stage, levels = unique(stage)), expr, group = age, fill = factor(age, levels = age_names))) + 
                        ggplot2::geom_line(ggplot2::aes(color = factor(age, levels = age_names)), size = 3) +
                        ggplot2::labs(x = xlab, y = ylab, title = main, colour = legendName) +
                        ggplot2::theme_minimal() +
                        ggplot2::theme(
                                title            = ggplot2::element_text(size = 18, face = "bold"),
                                legend.title     = ggplot2::element_text(size = 14, face = "bold"),
                                legend.text      = ggplot2::element_text(size = 18, face = "bold"),
                                axis.title       = ggplot2::element_text(size = 18, face = "bold"),
                                axis.text.y      = ggplot2::element_text(size = 18, face = "bold"),
                                axis.text.x      = ggplot2::element_text(size = 18, face = "bold"),
                                panel.background = ggplot2::element_blank(),
                                strip.text.x     = ggplot2::element_text(
                                        size           = 18,
                                        colour         = "black",
                                        face           = "bold"
                                )
                        ) +
                        ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = y.ticks)) + 
                        ggplot2::scale_colour_manual(values = custom.myTAI.cols(nrow(mMatrix))) 
                
                
                if (!is.null(modules)) {
                        p <- p + ggplot2::geom_rect(data = mMatrix,ggplot2::aes(
                                xmin = modules[[2]][1],
                                xmax = modules[[2]][length(modules[[2]])],
                                ymin = min(MeanValsMatrix) - (min(MeanValsMatrix) / 50),
                                ymax = Inf), fill = "#4d004b", alpha = alpha)  
                }
                p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1))
                return(p)
        }
        
        if (length(Groups) == 2) {
                
                mMatrixGroup1 <- dplyr::filter(mMatrix, age %in% Groups[[1]])
                mMatrixGroup2 <- dplyr::filter(mMatrix, age %in% Groups[[2]])
                
                p1 <- ggplot2::ggplot(mMatrixGroup1, ggplot2::aes( factor(stage, levels = unique(stage)), expr, group = age, fill = factor(age, levels = age_names[Groups[[1]]]))) + 
                        ggplot2::geom_line(ggplot2::aes(color = factor(age, levels = age_names[Groups[[1]]])), size = 3) +
                        ggplot2::labs(x = xlab, y = ylab, title = main, colour = legendName) +
                        ggplot2::theme_minimal() +
                        ggplot2::theme(
                                title            = ggplot2::element_text(size = 18, face = "bold"),
                                legend.title     = ggplot2::element_text(size = 14, face = "bold"),
                                legend.text      = ggplot2::element_text(size = 18, face = "bold"),
                                axis.title       = ggplot2::element_text(size = 18, face = "bold"),
                                axis.text.y      = ggplot2::element_text(size = 18, face = "bold"),
                                axis.text.x      = ggplot2::element_text(size = 18, face = "bold"),
                                panel.background = ggplot2::element_blank(),
                                strip.text.x     = ggplot2::element_text(
                                        size           = 18,
                                        colour         = "black",
                                        face           = "bold"
                                )
                        ) +
                        ggplot2::scale_colour_manual(values = custom.myTAI.cols(nrow(mMatrix))[Groups[[1]]]) +
                        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = -90, hjust = 0))
                if (!adjust.range) {
                        
                        p1 <- p1 + ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = y.ticks))
                }
                
                if (!is.null(modules)) {
                        p1 <- p1 + ggplot2::geom_rect(data = mMatrixGroup1, ggplot2::aes(
                                xmin = modules[[2]][1],
                                xmax = modules[[2]][length(modules[[2]])],
                                ymin = min(MeanValsMatrix),
                                ymax = Inf), fill = "#4d004b", alpha = alpha)  
                }
                
                p2 <- ggplot2::ggplot(mMatrixGroup2, ggplot2::aes( factor(stage, levels = unique(stage)), expr, group = age, fill = factor(age, levels = age_names[Groups[[2]]]))) + 
                        ggplot2::geom_line(ggplot2::aes(color = factor(age, levels = age_names[Groups[[2]]])), size = 3) +
                        ggplot2::labs(x = xlab, y = ylab, title = main, colour = legendName) +
                        ggplot2::theme_minimal() +
                        ggplot2::theme(
                                title            = ggplot2::element_text(size = 18, face = "bold"),
                                legend.title     = ggplot2::element_text(size = 14, face = "bold"),
                                legend.text      = ggplot2::element_text(size = 18, face = "bold"),
                                axis.title       = ggplot2::element_text(size = 18, face = "bold"),
                                axis.text.y      = ggplot2::element_text(size = 18, face = "bold"),
                                axis.text.x      = ggplot2::element_text(size = 18, face = "bold"),
                                panel.background = ggplot2::element_blank(),
                                strip.text.x     = ggplot2::element_text(
                                        size           = 18,
                                        colour         = "black",
                                        face           = "bold"
                                )
                        ) + 
                        ggplot2::scale_colour_manual(values = custom.myTAI.cols(nrow(mMatrix))[Groups[[2]]]) +
                        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = -90, hjust = 0))
                
                if (!adjust.range) {
                        
                        p2 <- p2 + ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = y.ticks)) 
                }
                
                
                if (adjust.range){
                        p1 <- p1 + ggplot2::scale_y_continuous(limits = c(min(MeanValsMatrix), max(MeanValsMatrix)), breaks = scales::pretty_breaks(n = y.ticks))
                        
                        p2 <- p2 + ggplot2::scale_y_continuous(limits = c(min(MeanValsMatrix), max(MeanValsMatrix)), breaks = scales::pretty_breaks(n = y.ticks))    
                }
                
                if (!is.null(modules)) {
                        p2 <- p2 + ggplot2::geom_rect(data = mMatrixGroup2,ggplot2::aes(
                                xmin = modules[[2]][1],
                                xmax = modules[[2]][length(modules[[2]])],
                                ymin = min(MeanValsMatrix),
                                ymax = Inf), fill = "#4d004b", alpha = alpha)  
                }
                p1 <- p1 + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1))
                p2 <- p2 + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1))
                return(gridExtra::grid.arrange(p1, p2, ncol = 2))
        }
}
