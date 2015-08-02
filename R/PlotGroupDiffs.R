#' @title Plot the significant differences between gene expression distributions of PS or DS groups
#' @description This function performs a statistical test to quantify the statistical significance between
#' the global expression level distributions of groups of PS or DS. In therefore, allows users to investigate
#' significant groups of PS or DS that significantly differ in their gene expression level distibution
#' within specific developmental stages or experiments.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param Groups a list containing the phylostrata or divergence strata that correspond 
#' to the same phylostratum class or divergence class.
#' For ex. evolutionary old phylostrata: PS1-3 (Class 1) 
#' and evolutionary young phylostrata: PS4-12 (Class 2). In this case, 
#' the list could be assigned as, \code{Groups} = list(c(1:3), c(4:12)).
#' @param legendName a character string specifying whether "PS" or "DS" are used to compute relative expression profiles. 
#' @param stat.test the statistical test to quantify PS or DS group differences.
#' @param plot.p.vals a logical value indicating whether the plot should be drawn or only the p-value should be returned without drawing the P-Value plot.
#' @param ... additional plot parameters.
#' @author Hajk-Georg Drost
#' @details 
#' The purpose of this function is to detect groups of PS or DS that significantly differ in their gene expression
#' level distributions on a global (transcriptome) level. Since relative expression levels (\code{\link{PlotRE}}) or
#' PS or DS specific mean expression levels (\code{\link{PlotMeans}}) are biased by highly expressed genes,
#' this function allows users to objectively test the significant difference of transcriptome expression between
#' groups of PS or DS in a specific developmental stage or experiment.
#' @examples 
#' 
#' data(PhyloExpressionSetExample)
#' 
#' PlotGroupDiffs(ExpressionSet = PhyloExpressionSetExample,
#'                Groups        = list(group_1 = 1:3,group_2 = 4:12),
#'                legendName    = "PS",
#'                type          = "b",
#'                lwd           = 6,
#'                xlab          = "Ontogeny")
#'                
#'                
#' # only receive the p-values without the corresponding plot               
#' PlotGroupDiffs(ExpressionSet = PhyloExpressionSetExample,
#'                Groups        = list(group_1 = 1:3,group_2 = 4:12),
#'                legendName    = "PS",
#'                plot.p.vals   = FALSE,
#'                type          = "b",
#'                lwd           = 6,
#'                xlab          = "Ontogeny")
#' 
#' @seealso \code{\link{PlotMeans}}, \code{\link{PlotRE}}, \code{\link{PlotBarRE}}, \code{\link{PlotCategoryExpr}}
#' @export

PlotGroupDiffs <- function(ExpressionSet,
                           Groups      = NULL,
                           legendName  = NULL,
                           stat.test   = "wilcox.test",
                           plot.p.vals = TRUE, ...){
        
        is.ExpressionSet(ExpressionSet)
       
        if (is.null(Groups))
                stop ("Your Groups list does not store any items.")
        
        if (is.null(legendName))
                stop ("Please specify the type of ExpressionSet you are working with: legendName = 'PS' or 'DS'.")
        
        ### getting the PS names available in the given expression set
        age_names <- as.character(names(table(ExpressionSet[ , 1])))
        ncols <- ncol(ExpressionSet)
        nStages <- ncols - 2
        
        # test whether all group elements are available in the age vector
        if (!all(unlist(Groups) %in% as.numeric(age_names)))
                stop ("There are items in your Group elements that are not available in the age column of your ExpressionSet.") 
        
        if (!is.element(stat.test, c("wilcox.test")))
                stop ("Unfortunately the statistical test '",stat.test,"' is not implemented in this function.")
        
        if (stat.test == "wilcox.test"){
                
                if (length(Groups) > 2)
                        stop ("To perform a pairwise wilcox.test you can only specify two sets of group elements.")
                
                p.val.stages <- vector("numeric", length = nStages)
                
                for (i in 1:nStages){
                        
                        p.val.stages[i] <- wilcox.test(ExpressionSet[which(ExpressionSet[ , 1] %in% Groups[[1]]), i + 2],
                                                       ExpressionSet[which(ExpressionSet[ , 1] %in% Groups[[2]]), i + 2])$p.value 
                }
        }
        
        if (plot.p.vals){
                # define arguments for different graphics functions
                plot.args <- c("lwd","col","lty","xlab","cex.lab","main","type")
                axis.args <- c("las", "cex.axis")
                legend.args <- c("border","angle","density","box.lwd","cex")
                dots <- list(...)
                ellipsis.names <- names(dots)
                
                
                do.call(graphics::plot,c(list(x = p.val.stages, xaxt = "n", ylab = "P-Value"),dots[!is.element(names(dots),c(axis.args,legend.args))]))
                        do.call(graphics::axis,c(list(side = 1,at = seq(1,nStages,1), labels = names(ExpressionSet)[3:ncols]), 
                                                 dots[!is.element(names(dots),c(plot.args,legend.args))]))   
                        
                 # do.call(graphics::legend,c(list(x = "top",legend = c(paste0(legendName,Groups[[1]]," "),paste0(legendName,Groups[[2]]," ")), bty = "n", ncol = 2, col = c("black","blue")),dots[!is.element(names(dots),c(axis.args,plot.args))]))
                 
        }
        
        return(p.val.stages)
}






