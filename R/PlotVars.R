#' @title Plot Variance of Expression Profiles
#' @description This function computes for each age category the corresponding variance expression profile.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param Groups a list containing the age categories for which variance expression levels shall be drawn.
#' For ex. evolutionary users can compare old phylostrata: PS1-3 (Class 1) and evolutionary young phylostrata: PS4-12 (Class 2). 
#' In this example, the list could be assigned as, \code{Groups = list(c(1:3), c(4:12))}. 
#' The group options is limited to 2 Groups.
#' @param legendName a character string specifying the legend title.
#' @param xlab label of x-axis.
#' @param ylab label of y-axis.
#' @param main main text.
#' @param y.ticks number of ticks that shall be drawn on the y-axis.
#' @param adjust.range logical indicating whether or not the y-axis scale shall be adjusted to the same range in case two groups are specified. Default is \code{adjust.range = TRUE}.
#' @details 
#' 
#' This plot may be useful to compare the absolute variance expression        
#' levels of each age category across stages.
#'
#'In different developmental processes different phylostratum or divergence-stratum
#' classes might be more expressed than others, hence contributing more to the overall
#' phylotranscriptomics pattern (\code{\link{TAI}} or \code{\link{TDI}}).
#' This plot can help to identify the phylostratum or divergence-stratum classes 
#' that contributes most to the overall transcriptome of the given developmental process.
#' @return a plot showing variance expression profiles of each age category.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{PlotBarRE}}, \code{\link{RE}}, \code{\link{REMatrix}}, \code{\link{PlotRE}}
#' @examples
#' ### Example using a PhyloExpressionSet
#' ### and DivergenceExpressionSet
#' # load PhyloExpressionSet
#' data(PhyloExpressionSetExample)
#'
#' # load PhyloExpressionSet
#' data(DivergenceExpressionSetExample)
#'
#' # plot evolutionary old PS (PS1-3) vs evolutionary young PS (PS4-12)
#' PlotVars(PhyloExpressionSetExample,
#'           Groups = list(c(1:3), c(4:12)), 
#'           legendName = "PS",
#'           adjust.range = TRUE)
#'
#' # if users wish to not adjust the y-axis scale when 
#' # 2 groups are selected they can specify: adjust.range = FALSE
#' PlotVars(PhyloExpressionSetExample,
#'           Groups = list(c(1:3), c(4:12)), 
#'           legendName = "PS",
#'           adjust.range = FALSE)
#'           
#'           
#' # plot conserved DS (DS1-5) vs divergent DS (PS6-10)
#' # NOTE: DS are always defined in the range 1, 2, ... , 10.
#' # Hence, make sure that your groups are within this range!
#' PlotVars(DivergenceExpressionSetExample,
#'           Groups = list(c(1:5), c(6:10)), 
#'           legendName = "DS",
#'           adjust.range = TRUE)
#'
#' @export 

PlotVars <- function(ExpressionSet,
                      Groups     = NULL,
                      legendName = "age",
                      xlab = "Ontogeny",
                      ylab = "Variance(Expression Level)",
                      main = "",
                      y.ticks = 10,
                      adjust.range = TRUE) 
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
        
        MeanValsMatrix <- age.apply(ExpressionSet, function(x) apply(x, 2, stats::var))
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
                p1 <- p1 + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1))
                p2 <- p2 + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1))
                return(gridExtra::grid.arrange(p1, p2, ncol = 2))
        }
}


