#' @title Plot Cumuative Transcriptome Index
#' @description This function computes the cumulative contribution of each Phylostratum or Divergence Stratum to the global \code{\link{TAI}} or \code{\link{TDI}} profile.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param legendName a character string specifying whether "PS" or "DS" are used to compute relative expression profiles.
#' @param xlab label of x-axis.
#' @param ylab label of y-axis.
#' @param main main title.
#' @param y.ticks a numeric value specifying the number of ticks to be drawn on the y-axis.
#' @details
#' Introduced by Domazet-Loso and Tautz (2010), this function allows users to visualize the cumulative contribution of each Phylostratum or Divergence Stratum to the global Transcriptome Age Index or Transcriptome Divergence Index profile to quantify how each Phylostratum or Divergence Stratum influences the profile of the global TAI or TDI pattern. 
#' 
#' @author Hajk-Georg Drost
#' @examples
#' 
#'  data(PhyloExpressionSetExample)
#'  data(DivergenceExpressionSetExample)
#'  
#'  # visualize phylostratum contribution to global TAI
#'  PlotContribution(PhyloExpressionSetExample, legendName = "PS")
#'  
#'  # visualize divergence stratum contribution to global TDI
#'  PlotContribution(DivergenceExpressionSetExample, legendName = "DS")
#'  
#' @references
#' 
#' Domazet-Loso T. and Tautz D. (2010). A phylogenetically based transcriptome age index mirrors ontogenetic divergence patterns. Nature (468): 815-818.
#'    
#' @seealso \code{\link{pTAI}}, \code{\link{pTDI}}, \code{\link{TAI}}, \code{\link{TDI}}, \code{\link{PlotSignature}}
#' @export         

PlotContribution <- function(ExpressionSet, 
                             legendName  = NULL,
                             xlab = "Ontogeny",
                             ylab = "Transcriptome Index",
                             main = "",
                             y.ticks     = 10){
        
        if(is.null(legendName))
                stop("Please specify whether your input ExpressionSet stores 'PS' or 'DS'.", call. = FALSE)
    
        ExpressionSet <- as.data.frame(ExpressionSet)
        
        DS <- par_value <-  PS <- stage <- NULL
        
        is.ExpressionSet(ExpressionSet)
        
        ncols <- ncol(ExpressionSet)
        nPS <- length(table(ExpressionSet[ , 1]))
        
        # define contribution matrix
        contrMatrix <- matrix(NA_real_,ncol = ncols,nrow = nPS)
        contrMatrix <- pTAI(ExpressionSet)

        if (legendName == "PS") {
                
                contrMatrix <- data.frame(PS = paste0("PS",1:nrow(contrMatrix)),contrMatrix, 
                                          stringsAsFactors = FALSE)
                contrMatrix <- tibble::as_tibble(contrMatrix)
                contrMatrix <- reshape2::melt(contrMatrix, id.vars = "PS")
                colnames(contrMatrix)[2:3] <- c("stage", "par_value")
                
                p <- ggplot2::ggplot(contrMatrix, ggplot2::aes( factor(stage, levels = unique(stage)), par_value, group = PS, fill = factor(PS, levels = paste0("PS",1:nrow(contrMatrix))))) + 
                        ggplot2::geom_line(ggplot2::aes(color = factor(PS, levels = paste0("PS",1:nrow(contrMatrix)))), size = 3) +
                        ggplot2::labs(x = xlab, y = ylab, title = main, colour = "PS") +
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
                        ggplot2::scale_colour_manual(values = custom.myTAI.cols(nrow(contrMatrix)))
        }
        
        
        if (legendName == "DS") {
                
                contrMatrix <- data.frame(DS = paste0("DS",1:nrow(contrMatrix)),contrMatrix, 
                                          stringsAsFactors = FALSE)
                contrMatrix <- tibble::as_tibble(contrMatrix)
                contrMatrix <- reshape2::melt(contrMatrix, id.vars = "DS")
                colnames(contrMatrix)[2:3] <- c("stage", "par_value")
                
                p <- ggplot2::ggplot(contrMatrix, ggplot2::aes( factor(stage, levels = unique(stage)), par_value, group = DS, fill = factor(DS, levels = paste0("DS",1:nrow(contrMatrix))))) + 
                        ggplot2::geom_line(ggplot2::aes(color = factor(DS, levels = paste0("DS",1:nrow(contrMatrix)))), size = 3) +
                        ggplot2::labs(x = xlab, y = ylab, title = main, colour = "DS") +
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
                        ggplot2::scale_colour_manual(values = custom.myTAI.cols(nrow(contrMatrix)))
        }
        
        p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 1,hjust = 1))
        return(p)
}







