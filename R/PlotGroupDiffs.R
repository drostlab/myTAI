#' @title Plot the significant differences between gene expression distributions of PS or DS groups
#' @description This function performs a test to quantify the statistical significance between
#' the global expression level distributions of groups of PS or DS. It therefore, allows users to investigate
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
#' @param col colors for the two box plots representing the expression level distributions of selected PS/DS groups.
#' @param plot.type the type of plot that shall be drawn to visualized the difference in PS/DS group specific expression .
#' @param gene.set a character vector storing the gene ids for which group specific differences shall be statistically quantified.
#' @param ... additional plot parameters.
#' @author Hajk-Georg Drost
#' @details 
#' The purpose of this function is to detect groups of PS or DS that significantly differ in their gene expression
#' level distributions on a global (transcriptome) level. Since relative expression levels (\code{\link{PlotRE}}) or
#' PS or DS specific mean expression levels (\code{\link{PlotMeans}}) are biased by highly expressed genes,
#' this function allows users to objectively test the significant difference of transcriptome expression between
#' groups of PS or DS in a specific developmental stage or experiment.
#' 
#' In particular, this function divides (for each developmental stage separately) the gene expression levels into two groups: Group1 = genes deriving from selected PS/DS in group 1 and 
#' Group2 = genes deriving from selected PS/DS in group 2. Within each stage the expression level distributions between group 1 and group 2 are statistically quantified using a \code{\link{wilcox.test}}.
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
#'                
#' # quantify the significant difference of a selected set of genes
#' # only receive the p-values without the corresponding plot
#' set.seed(123)
#' ExampleGeneSet <- sample(PhyloExpressionSetExample[ , 2],5000)
#'                               
#' PlotGroupDiffs(ExpressionSet = PhyloExpressionSetExample,
#'                Groups        = list(group_1 = 1:3,group_2 = 4:12),
#'                legendName    = "PS",
#'                plot.p.vals   = FALSE,
#'                gene.set      = ExampleGeneSet)                 
#' 
#' 
#' # plot differences as boxplot for each developmental stage
#' PlotGroupDiffs(ExpressionSet = tf(PhyloExpressionSetExample,log2),
#'                Groups        = list(group_1 = 1:3,group_2 = 4:12),
#'                legendName    = "PS",
#'                plot.type     = "boxplot")
#' 
#' 
#' @seealso \code{\link{PlotMeans}}, \code{\link{PlotRE}}, \code{\link{PlotBarRE}}, \code{\link{PlotCategoryExpr}}, \code{\link{GroupDiffs}}
#' @export

PlotGroupDiffs <- function(ExpressionSet,
                           Groups      = NULL,
                           legendName  = NULL,
                           stat.test   = "wilcox.test",
                           col         = c("turquoise3","magenta3"),
                           plot.type   = NULL,
                           gene.set    = NULL, ...){
        
        is.ExpressionSet(ExpressionSet)
       
        if (is.null(Groups))
                stop ("Your Groups list does not store any items.")
        
        if (is.null(legendName))
                stop ("Please specify the type of ExpressionSet you are working with: legendName = 'PS' or 'DS'.")
        if (length(col) > 2)
                stop ("Please enter only two colors for the two groups.")
        
        if (!is.element(stat.test,c("wilcox.test")))
                stop (stat.test, " is not implemented in this function.")
       
        if (!is.null(plot.type)) 
                if (!is.element(plot.type,c("p-vals","boxplot")))
                        stop ("Please select a plot.type that is supported by this function.")
        
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
                
        if (!is.null(gene.set)){
                
                ExpressionSet <- SelectGeneSet(ExpressionSet = ExpressionSet, gene.set = gene.set)
        }        
                
                p.val.stages <- vector("numeric", length = nStages)
                GroupCategoryList <- 1
                group.names <- c("Group1","Group2")  
                
                if (!is.null(plot.type))
                        if (plot.type == "boxplot")
                                graphics::par(mfrow = grDevices::n2mfrow(nStages))
                
                for (i in 1:nStages){
                        
                        Group1 <- ExpressionSet[which(ExpressionSet[ , 1] %in% Groups[[1]]), i + 2]
                        Group2 <- ExpressionSet[which(ExpressionSet[ , 1] %in% Groups[[2]]), i + 2]
                        
                        p.val.stages[i] <- stats::wilcox.test(ExpressionSet[which(ExpressionSet[ , 1] %in% Groups[[1]]), i + 2],
                                                       ExpressionSet[which(ExpressionSet[ , 1] %in% Groups[[2]]), i + 2])$p.value 
                        if (!is.null(plot.type)){
                                if (plot.type == "boxplot"){
                                        dataList <- lapply(group.names, get, envir=environment())
                                        names(dataList) <- group.names
                                        graphics::boxplot(dataList,xlab = "Groups", ylab = "Expression Level", main = paste0(names(ExpressionSet)[i + 2],"  ( P = ",format(p.val.stages[i],digits = 2)," )"), col = col)  
                                }  
                        }
                }
                
                p.val.stages <- t(as.data.frame(p.val.stages))
                colnames(p.val.stages) <- names(ExpressionSet)[3:ncols]
                rownames(p.val.stages) <- paste0("p.value ( ",stat.test," )")
        }
        
        
        if (!is.null(plot.type)){
                if (plot.type == "p-vals"){
                        
                        # define arguments for different graphics functions
                        plot.args <- c("lwd","col","lty","xlab","cex.lab","main","type")
                        axis.args <- c("las", "cex.axis")
                        legend.args <- c("border","angle","density","box.lwd","cex")
                        dots <- list(...)
                        ellipsis.names <- names(dots)
                        
                        
                        do.call(graphics::plot,c(list(x = 1:nStages,y = p.val.stages, xaxt = "n", ylab = "P-Value"),dots[!is.element(names(dots),c(axis.args,legend.args))]))
                        do.call(graphics::axis,c(list(side = 1,at = seq(1,nStages,1), labels = names(ExpressionSet)[3:ncols]), 
                                                 dots[!is.element(names(dots),c(plot.args,legend.args))]))   
                        
                        # do.call(graphics::legend,c(list(x = "top",legend = c(paste0(legendName,Groups[[1]]," "),paste0(legendName,Groups[[2]]," ")), bty = "n", ncol = 2, col = c("black","blue")),dots[!is.element(names(dots),c(axis.args,plot.args))]))
                        
                }
        }
        
        return(p.val.stages)
}






