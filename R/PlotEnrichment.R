#' @title Plot the Phylostratum or Divergence Stratum Enrichment of a given Gene Set
#' @description This function computes and visualizes the significance of enriched (over or underrepresented) Phylostrata or Divergence Strata within an input \code{test.set}.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object (in case \code{only.map = FALSE}).
#' @param test.set a character vector storing the gene ids for which PS/DS enrichment analyses should be performed.
#' @param use.only.map a logical value indicating whether instead of a standard \code{ExpressionSet} only a \code{Phylostratigraphic Map} or \code{Divergene Map} is passed to this function.
#' @param measure a character string specifying the measure that should be used to quantify over and under representation of PS/DS. Measures can either be \code{measure = "foldchange"} (odds) or \code{measure = "log-foldchange"} (log-odds).
#' @param complete.bg a logical value indicating whether the entire background set
#'  of the input \code{ExpressionSet} should be considered when performing Fisher's exact
#'  test (\code{complete.bg = TRUE}) or whether genes that are stored in \code{test.set}
#'  should be excluded from the background set before performing Fisher's exact test
#'  (\code{complete.bg = FALSE}).
#' @param legendName a character string specifying whether "PS" or "DS" are used to compute relative expression profiles.
#' @param over.col color of the overrepresentation bars.
#' @param under.col color of the underrepresentation bars.
#' @param epsilon a small value to shift values by epsilon to avoid log(0) computations.
#' @param plot.bars a logical value specifying whether or not bars should be visualized or whether only \code{p.values} and \code{enrichment.matrix} should be returned.
#' @param cex.legend the \code{cex} value for the legend.
#' @param cex.asterisk the \code{cex} value for the asterisk.
#' @param p.adjust.method correction method to adjust p-values for multiple comparisons (see \code{\link{p.adjust}} for possible methods). E.g., \code{p.adjust.method = "BH"} (Benjamini & Hochberg (1995)) or \code{p.adjust.method = "bonferroni"} (Bonferroni correction).
#' @param ... default graphics parameters.
#' @details 
#' 
#' This \emph{Phylostratum} or \emph{Divergence Stratum} enrichment analysis is motivated 
#' by Sestak and Domazet-Loso (2015) who perform \emph{Phylostratum} or \emph{Divergence Stratum} enrichment analyses to correlate organ evolution with the origin of organ specific genes. 
#' 
#' In detail this function takes the \emph{Phylostratum} or \emph{Divergence Stratum} distribution of all genes stored in the input \code{ExpressionSet} as background set and
#' the \emph{Phylostratum} or \emph{Divergence Stratum} distribution of the \code{test.set} and performes a \code{\link{fisher.test}} for each \emph{Phylostratum} or \emph{Divergence Stratum} to quantify the statistical significance of over- or underrepresentated \emph{Phylostrata} or \emph{Divergence Strata} within the set of selected \code{test.set} genes. 
#' 
#' To visualize the odds or log-odds of over or underrepresented genes within the \code{test.set} the following procedure is performed:
#' 
#' \itemize{
#' \item N_ij denotes the number of genes in group j and deriving from PS i, with \emph{i = 1, .. , n} and where \emph{j = 1} denotes the background set and \emph{j = 2} denotes the \code{test.set}
#' \item N_i. denotes the total number of genes within PS i
#' \item N_.j denotes the total number of genes within group j
#' \item N_.. is the total number of genes within all groups j and all PS i
#' \item f_ij = N_ij / N_.. and g_ij = f_ij / f_.j denote relative frequencies between groups
#' \item f_i. denotes the between group sum of f_ij
#' }
#' 
#'  The result is the fold-change value (odds) denoted as C = g_i2 / f_i. which is visualized above and below zero. 
#'  
#'  In case a large number of Phylostrata or Divergence Strata is included in the input 
#'  \code{ExpressionSet}, p-values returned by \code{PlotEnrichment} should be adjusted for
#'  multiple comparisons which can be done by specifying the \code{p.adjust.method} argument.
#'  
#' @author Hajk-Georg Drost
#' @references
#' 
#' Sestak and Domazet-Loso (2015). Phylostratigraphic Profiles in Zebrafish Uncover Chordate Origins of the Vertebrate Brain. Mol. Biol. Evol. 32(2): 299-312. 
#' 
#' @examples
#' 
#' data(PhyloExpressionSetExample)
#' 
#' set.seed(123)
#' test_set <- sample(PhyloExpressionSetExample[ , 2],10000)
#' 
#' ## Examples with complete.bg = TRUE
#' ## Hence: the entire background set of the input ExpressionSet is considered 
#' ## when performing Fisher's exact test 
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
#'                
#' ## Examples with complete.bg = FALSE
#' ## Hence: the test.set genes are excluded from the background set before
#' ## Fisher's exact test is performed
#'      
#'                                        
#' # measure: log-foldchange
#' PlotEnrichment(ExpressionSet = PhyloExpressionSetExample,
#'                test.set      = test_set ,
#'                 complete.bg  = FALSE,
#'                legendName    = "PS", 
#'                measure       = "log-foldchange")
#'                
#'                
#' # measure: foldchange
#' PlotEnrichment(ExpressionSet = PhyloExpressionSetExample,
#'                test.set      = test_set , 
#'                complete.bg   = FALSE,
#'                legendName    = "PS", 
#'                measure       = "foldchange")     
#'                
#' @seealso \code{\link{EnrichmentTest}}, \code{\link{fisher.test}}         
#' @export

PlotEnrichment <- function(ExpressionSet,
                           test.set,
                           use.only.map    = FALSE,
                           measure         = "log-foldchange",
                           complete.bg     = TRUE,
                           legendName      = "",
                           over.col        = "steelblue",
                           under.col       = "midnightblue",
                           epsilon         = 0.00001,
                           cex.legend      = 1,
                           cex.asterisk    = 1,
                           plot.bars       = TRUE,
                           p.adjust.method = NULL, ...){
        
        ExpressionSet <- as.data.frame(ExpressionSet)
        
        # check for correct data input (ExpressionSet or Phylomap/Divergencemap) 
        if(!use.only.map){
                
                if (ncol(ExpressionSet) < 3)
                        warning("Are you sure that you are using an ExpressionSet and not a Phylomap/Divergencemap ?")
                
                is.ExpressionSet(ExpressionSet)
        } else {
                if(!is.numeric(ExpressionSet[ , 1]) | !is.character(ExpressionSet[ , 2]))
                        stop("Please pass a standard Phylo Map or Divergence Map to this function...")
        }
                
        if (!is.element(measure , c("log-foldchange","foldchange")))
                stop("Please select a measure which is supported by this function...")
       
        if (length(test.set) > nrow(ExpressionSet))
                stop("Your input GeneID vector stores more elements than are available in your ExpressionSet object...")
        
        age.table <- table(ExpressionSet[ , 1])
        nPS <- length(age.table)
        
        MatchedGeneIDs <- stats::na.omit(match(tolower(unlist(test.set)),tolower(unlist(ExpressionSet[ , 2]))))
        
        if (length(MatchedGeneIDs) == 0)
                stop("None of your input gene ids could be found in the ExpressionSet.")
        
        age.distr.test.set <- ExpressionSet[MatchedGeneIDs , 1:2]
        
        # exclude test.set genes from the background set
        # before performing Fisher's exact test
        if (complete.bg)
                ExpressionSet <- ExpressionSet[-MatchedGeneIDs , ]

        if (length(age.distr.test.set[ , 2]) != length(test.set))
                warning ("Only ",length(age.distr.test.set[ , 2]), " out of your ",length(test.set)," gene ids could be found in the ExpressionSet.")
        
        FactorBackgroundSet <- factor(ExpressionSet[ , 1],levels = names(age.table))
        FactorTestSet <- factor(age.distr.test.set[ , 1],levels = names(age.table))
        
        # number of genes in group j having PS i
        N_ij <- cbind(table(FactorBackgroundSet),table(FactorTestSet));
        # naming the 2 groups : group1 = "Background" ; group2 = "Test"
        colnames(N_ij) <- c("Group 1: Background","Group 2: Test")
        
        enrichment.p_vals <- vector("numeric",nPS)
        enrichment.p_vals <- sapply(1:nPS,
                                    function(index) stats::fisher.test(get.contingency.tbl(N_ij,index))$p.value,
                                    simplify = "array")
        
        names(enrichment.p_vals) <- paste0(legendName,names(age.table))
        
        # number of all genes in N_ij
        N_dot_dot <- sum(N_ij)
        
        # rel. freq. over all elements in both groups
        f_ij <- N_ij / N_dot_dot
        
        # between group sum
        f_i_dot <- rowSums(f_ij)
        
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
                group_1 <- vector("numeric",nPS)
                group_2 <- vector("numeric",nPS)
                
                for(i in 1:nPS){
                        if(f_i_dot[i] > 0){
                                if(g_ij[i , 1] > 0){
                                        group_1[i] <- (log2(g_ij[i , 1] + epsilon) - log2(f_i_dot[i] + epsilon))
                                } else {
                                        group_1[i] <- 0
                                }
                                
                                if(g_ij[i , 2] > 0){
                                        group_2[i] <- (log2(g_ij[i , 2] + epsilon) - log2(f_i_dot[i] + epsilon))
                                } else {
                                        group_2[i] <- 0
                                }
                                
                        } else {
                                group_1[i] <- 0
                                group_2[i] <- 0
                        }
                } 
                ResultMatrix <- cbind(group_1, group_2)
                
                UpRegulated <- which(ResultMatrix[ , 2] >= 0)
                DownRegulated <- which(ResultMatrix[ , 2] < 0)
        }
        
        if (measure == "foldchange"){
                # epsilon is a shift operator to omit log2(0) = -Inf vals
                # ResultMatrix[,1] shows over or underrepresentation of group 1 by factor X compared to group 2
                # ResultMatrix[,2] shows over or underrepresentation of group 2 by factor X compared to group 1
                group_1 <- vector("numeric",nPS)
                group_2 <- vector("numeric",nPS)
                
                for(i in 1:nPS){
                        if(f_i_dot[i] > 0){
                                if(g_ij[i , 1] > 0){
                                        group_1[i] <- (g_ij[i , 1] / f_i_dot[i])
                                } else {
                                        group_1[i] <- 0
                                }
                                        
                                if(g_ij[i , 2] > 0){
                                        group_2[i] <- (g_ij[i , 2] / f_i_dot[i])
                                } else {
                                        group_2[i] <- 0
                                }
                                
                        } else {
                                group_1[i] <- 0
                                group_2[i] <- 0
                        }
                } 
                
                ResultMatrix <- cbind(group_1 ,group_2 )
                
                # detect up and down regulated age classes 
                ResultMatrix[which((ResultMatrix[ , 2] < 1) & (ResultMatrix[ , 2]!=0)), 2] <- (-(1 / (ResultMatrix[which((ResultMatrix[ , 2] < 1) & (ResultMatrix[ , 2] !=0 )), 2])))
                # define a value for INF values
                ResultMatrix[which(ResultMatrix[ , 2] == 0), 2] <- 0
                UpRegulated <- which(ResultMatrix[ , 2] >= 1)
                DownRegulated <- which(ResultMatrix[ , 2] < 1)
        }
        
        colnames(ResultMatrix) <- c("BG_Set","Test_Set")
        rownames(ResultMatrix) <- paste0(legendName,names(age.table))
                
        if (plot.bars){
                
                dots <- list(...)
                ellipsis.names <- names(dots)
                
                # plot the results and only show the values (gi2 - f_i_dot) 
                # denoting under/over-representation of group2 = TestSet when compared
                # with group 1 = BackgroundSet
                
                bar.colors <- vector("character",nPS)
                bar.colors[UpRegulated] <- over.col
                bar.colors[DownRegulated] <- under.col
                
                max.lim.val <- max(abs(c(floor(min(ResultMatrix[ , 2]) - 0.1), ceiling(max(ResultMatrix[ , 2]) + 0.3))))
                ylim.range <- range(-max.lim.val,max.lim.val)
                
                if((length(ellipsis.names[grep("ylab",ellipsis.names)]) > 0) || (length(ellipsis.names[grep("xlab",ellipsis.names)]) > 0)){
                        
                        if (legendName == "PS")
                                age_xlab <- "Phylostratum"
                        
                        if (legendName == "DS")
                                age_xlab <- "Divergence Stratum"
                        
                        if (measure == "foldchange")
                                measure_ylab <- "Fold Change"
                        
                        if (measure == "log-foldchange")
                                measure_ylab <- "Log odds"
                        
                        barPlotFoldChanges <- graphics::barplot(ResultMatrix[ , 2],
                                                                beside    = FALSE,
                                                                col       = bar.colors,
                                                                lwd       = 4,
                                                                space     = 1,
                                                                names.arg = paste0(legendName,names(age.table)),
                                                                ylim      = ylim.range,
                                                                border    = "white",
                                                                xlab      = age_xlab, 
                                                                ylab      = measure_ylab, ...)
                        
                } else {
                        
                        barPlotFoldChanges <- graphics::barplot(ResultMatrix[ , 2],
                                                                beside    = FALSE,
                                                                col       = bar.colors,
                                                                lwd       = 4,
                                                                space     = 1,
                                                                names.arg = paste0(legendName,names(age.table)),
                                                                ylim      = ylim.range,
                                                                border    = "white", ...)
                }
                
                graphics::abline(h = 0, col = "black", lty = 1,lwd = 4)
                graphics::legend("topright",
                                 legend = c("over-represented","under-represented"),
                                 fill   = c(over.col,under.col),
                                 bty    = "n",
                                 cex    = cex.legend)
                
                
                pValNames <- vector("character",nPS)
                pValNames <- rep("",nPS)
                
                # adding multiple testing correction option
                if(!is.null(p.adjust.method)){
                        enrichment.p_vals <- stats::p.adjust(enrichment.p_vals,method = p.adjust.method, n = nPS)
                }
                
                pValNames[which(enrichment.p_vals <= 0.05)] <- "*"
                pValNames[which(enrichment.p_vals <= 0.005)] <- "**"
                pValNames[which(enrichment.p_vals <= 0.0005)] <- "***"
                pValNames[which(is.na(pValNames))] <- ""
                
                graphics::text(x      = barPlotFoldChanges,
                               y      = (max(ResultMatrix[ , 2]) + (ylim.range[2] / 5)),
                               labels = pValNames,
                               cex    = cex.asterisk)
        }
        
        return(list(p.values = enrichment.p_vals, enrichment.matrix = ResultMatrix))
}        
        
        
