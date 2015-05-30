#' @title Plot the Phylostratum or Divergence Stratum Enrichment of a given Gene Set
#' @description This function computes and visualizes the significance of enriched (over or underrepresented) Phylostrata or Divergence Strata within an input \code{test.set}.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param test.set a character vector storing the gene ids for which PS/DS enrichment analyses should be performed.
#' @param measure a character string specifying the measure that should be used to quantify over and under representation of PS/DS. 
#' @param legendName a character string specifying whether "PS" or "DS" are used to compute relative expression profiles.
#' @param over.col color of the overrepresentation bars.
#' @param under.col color of the underrepresentation bars.
#' @param epsilon a small value to shift values by epsilon to avoid log(0) computations.
#' @param plot.bars a logical value specifying whether or not bars should be visualized or whether only \code{p.values} and \code{enrichment.matrix} should be returned.
#' @param cex.legend the \code{cex} value for the legend.
#' @param cex.asterisk the \code{cex} value for the asterisk.
#' @param ... default graphics parameters.
#' @details This
#' @author Hajk-Georg Drost
#' @examples
#' 
#' data(PhyloExpressionSetExample)
#' 
#' set.seed(123)
#' test_set <- sample(PhyloExpressionSetExample[ , 2],1000)
#' 
#' # measure: log-foldchange
#' PlotEnrichment(ExpressionSet = PhyloExpressionSetExample,
#'                test.set      = test_set , 
#'                legendName    = "PS", 
#'                measure       = "log-foldchange")
#'                
#'                
#' # measure: foldchange
#' PlotEnrichment(ExpressionSet = PhyloExpressionSetExample,
#'                test.set      = test_set , 
#'                legendName    = "PS", 
#'                measure       = "foldchange")
#'                
#' @export

PlotEnrichment <- function(ExpressionSet,
                           test.set,
                           measure      = "log-foldchange",
                           legendName   = NULL,
                           over.col     = "steelblue",
                           under.col    = "midnightblue",
                           epsilon      = 0.00001,
                           cex.legend   = 1,
                           cex.asterisk = 1,
                           plot.bars    = TRUE, ...){
        
        is.ExpressionSet(ExpressionSet)
        
        if (!is.element(measure , c("log-foldchange","foldchange")))
                stop("Please select a measure which is supported by this function...")
        
        
#         if (!is.character(test.set))
#                 stop("Your input GeneIDs should be character strings...")
        
        if (length(test.set) > nrow(ExpressionSet))
                stop("Your input GeneID vector stores more elements than are available in your ExpressionSet object...")
        
        MatchedGeneIDs <- na.omit(match(tolower(test.set),tolower(ExpressionSet[ , 2])))
        age.distr.test.set <- ExpressionSet[MatchedGeneIDs , 1:2]
        
        if (length(age.distr.test.set[ , 2]) != length(test.set))
                warning(length(test.set) - length(age.distr.test.set[ , 2]), " out of ",length(test.set)," genes could not be found within the ExpressionSet object.")
        
        age.table <- table(ExpressionSet[ , 1])
        nPS <- length(age.table)
        
        FactorBackgroundSet <- factor(ExpressionSet[ , 1],levels = names(age.table))
        FactorTestSet <- factor(age.distr.test.set[ , 1],levels = names(age.table))
        
        # number of genes in group j having PS i
        N_ij <- cbind(table(FactorBackgroundSet),table(FactorTestSet));
        # naming the 2 groups : group1 = "Background" ; group2 = "Test"
        colnames(N_ij) <- c("Group 1: Background","Group 2: Test")
        
         
        #return(fisher.test(N_ij[ , 1],N_ij[ , 2],simulate.p.value = TRUE))
        
        enrichment.p_vals <- vector("numeric",nPS)
        enrichment.p_vals <- sapply(1:nPS,
                                    function(index) fisher.test(get.contingency.tbl(N_ij,index))$p.value,
                                    simplify = "array")
        
        names(enrichment.p_vals) <- paste0(legendName,1:nPS)
        
        # number of all genes in N_ij
        N_dot_dot <- sum(N_ij)
        
        # rel. freq. over all elements in both groups
        f_ij <- N_ij / N_dot_dot
        
        # between group sum
        f_i_dot <- colSums(f_ij)
        
        # within group sum
        f_dot_j <- colSums(f_ij)
        
        # defining the relative freq. Matrix of group 1 and group 2
        g_ij <- matrix(NA_real_,nrow(N_ij),2)
        g_ij[ ,1] <- N_ij[ , 1] / sum(N_ij[ , 1])
        g_ij[ ,2] <- N_ij[ , 2] / sum(N_ij[ , 2])
        colnames(g_ij) <- c("G1","G2")
        
        if (measure == "log-foldchange"){
                # epsilon is a shift operator to omit log2(0) = -Inf vals
                # ResultMatrix[,1] shows over or underrepresentation of group 1 by factor X compared to group 2
                # ResultMatrix[,2] shows over or underrepresentation of group 2 by factor X compared to group 1
                ResultMatrix <- cbind( (log2(g_ij[ , 1] + epsilon) - log2(f_i_dot + epsilon)), (log2(g_ij[ , 2] + epsilon) - log2(f_i_dot + epsilon)) )
                
                UpRegulated <- which(ResultMatrix[ , 2] >= 0)
                DownRegulated <- which(ResultMatrix[ , 2] < 0)
        }
        
        if (measure == "foldchange"){
                # epsilon is a shift operator to omit log2(0) = -Inf vals
                # ResultMatrix[,1] shows over or underrepresentation of group 1 by factor X compared to group 2
                # ResultMatrix[,2] shows over or underrepresentation of group 2 by factor X compared to group 1
                ResultMatrix <- cbind( (g_ij[ , 1] / f_i_dot),(g_ij[ , 2] / f_i_dot) )
                
                # detect up and down regulated age classes 
                ResultMatrix[which((ResultMatrix[ , 2] < 1) & (ResultMatrix[ , 2]!=0)), 2] <- (-(1 / (ResultMatrix[which((ResultMatrix[ , 2] < 1) & (ResultMatrix[ , 2] !=0 )), 2])))
                # define a value for INF values
                ResultMatrix[which(ResultMatrix[ , 2] == 0), 2] <- 0
                UpRegulated <- which(ResultMatrix[ , 2] >= 1)
                DownRegulated <- which(ResultMatrix[ , 2] < 1)
        }
        
        colnames(ResultMatrix) <- c("BG_Set","Test_Set")
        rownames(ResultMatrix) <- paste0(legendName,1:nPS)
        
        # plot the results and only show the values (gi2 - f_i_dot) 
        # denoting under/over-representation of group2 = TestSet when compared
        # with group 1 = BackgroundSet
        
        bar.colors <- vector("character",nPS)
        bar.colors[UpRegulated] <- over.col
        bar.colors[DownRegulated] <- under.col
        
        max.lim.val <- max(abs(c(floor(min(ResultMatrix[ , 2]) - 0.1), ceiling(max(ResultMatrix[ , 2]) + 0.3))))
        ylim.range <- range(-max.lim.val,max.lim.val)
        
        if(plot.bars){
                
                barPlotFoldChanges <- barplot(ResultMatrix[ , 2],
                                              beside    = FALSE,
                                              col       = bar.colors,
                                              lwd       = 4,
                                              space     = 1,
                                              names.arg = paste0(legendName,1:nPS),
                                              ylim      = ylim.range,
                                              border    = "white", ...)
                
                abline(h = 0,col = "black",lty = 1,lwd = 4)
                legend("topright",
                       legend = c("over-represented","under-represented"),
                       fill   = c(over.col,under.col),
                       bty    = "n",
                       cex    = cex.legend)
                
                
                pValNames <- vector("character",nPS)
                pValNames <- rep("",nPS)
                pValNames[which(enrichment.p_vals <= 0.05)] <- "*"
                pValNames[which(enrichment.p_vals <= 0.005)] <- "**"
                pValNames[which(enrichment.p_vals <= 0.0005)] <- "***"
                
                text(x      = barPlotFoldChanges,
                     y      = (max(ResultMatrix[ , 2]) + (ylim.range[2] / 5)),
                     labels = pValNames,
                     cex    = cex.asterisk)
        }
        
        return(list(p.values = enrichment.p_vals, enrichment.matrix = ResultMatrix))
        
}        
        
        
