#' @title Plot the Correlation Between Phylostrata and Divergence Strata
#' @description
#' This function plots the correlation coefficient between phylostratum values 
#' and divergence-stratum values of a given PhyloExpressionSet and DivergenceExpressionSet.
#'       
#' This function can be used to test whether a given PS distribution and DS distribution are 
#' linear correlated so that the independence of PS and DS can be assumed for 
#' subsequent analyses (Quint et al., 2012).
#' @param PhyloExpressionSet a standard PhyloExpressionSet object.
#' @param DivergenceExpressionSet a standard DivergenceExpressionSet object.
#' @param method a character string specifying the correlation method to cbe used, e.g. "pearson", "kendall", "spearman".  
#' @param linearModel a boolean value specifying whether a linear model should be
#' fitted to the data and furthermore, should be visualized in the corresponding plot.
#' @return a jitter-correlation-plot of PS and DS correlation.
#' @references  Quint M et al. (2012). \emph{A transcriptomic hourglass in plant embryogenesis}. Nature (490): 98-101.
#' Drost HG et al. (2015) \emph{Evidence for Active Maintenance of Phylotranscriptomic Hourglass Patterns in Animal and Plant Embryogenesis}. Mol Biol Evol. 32 (5): 1221-1231 doi:10.1093/molbev/msv012.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{cor}}
#' @examples 
#' 
#'  # read standard phylotranscriptomics data
#'  data(PhyloExpressionSetExample)
#'  data(DivergenceExpressionSetExample)
#'  
#' # plot the PS and DS correlation
#' PlotCorrelation(PhyloExpressionSetExample, DivergenceExpressionSetExample, 
#'                 method      = "pearson", 
#'                 linearModel = TRUE)
#' 
#' 
#' 
#' @export

PlotCorrelation <- function(PhyloExpressionSet,DivergenceExpressionSet,
                            method      = "pearson",
                            linearModel = FALSE)
{
        
        is.ExpressionSet(PhyloExpressionSet)
        is.ExpressionSet(DivergenceExpressionSet)
        
        if(!is.element(method, c("pearson", "kendall", "spearman")))
                stop("Please choose a correlation method that is supported by this function.")
        
        colnames(PhyloExpressionSet)[2] <- "GeneID"
        colnames(DivergenceExpressionSet)[2] <- "GeneID"
        nPS <- length(names(table(PhyloExpressionSet[ , 1])))
        nDS <- length(names(table(DivergenceExpressionSet[ , 1])))
        
        # convert ids to lower case
        PhyloExpressionSet[ , "GeneID"] <- tolower(PhyloExpressionSet[ , "GeneID"])
        DivergenceExpressionSet[ , "GeneID"] <- tolower(DivergenceExpressionSet[ , "GeneID"])
        
        PS_DS.Subset <- merge(PhyloExpressionSet[ , 1:2], DivergenceExpressionSet[ , 1:2],by = "GeneID")
        
        CorrelationCoefficient <- stats::cor(PS_DS.Subset[ , 2],PS_DS.Subset[ , 3],method = method)
        CorrCoeffasCharacter <- as.character(round(CorrelationCoefficient,3))
        
        nrows <- dim(PS_DS.Subset)[1]
        PS <- vector(mode = "numeric", length = nrows)
        DS <- vector(mode = "numeric", length = nrows)
        
        PS <- jitter(PS_DS.Subset[ , 2],1.5)
        DS <- jitter(PS_DS.Subset[ , 3],1.5)
        
        CorDF <- data.frame(PS = PS, DS = DS)
        
        if (!linearModel){
                
                res.plot <- ggplot2::ggplot(CorDF, ggplot2::aes(x=PS, y=DS)) + ggplot2::ggtitle(paste0(method," = ",CorrCoeffasCharacter)) +
                        ggplot2::geom_point(shape=1) +  ggplot2::theme_minimal()
                
        }
        
        if (linearModel){
         
                res.plot <- ggplot2::ggplot(CorDF, ggplot2::aes(x=PS, y=DS)) + ggplot2::ggtitle(paste0(method," = ",CorrCoeffasCharacter)) +
                        ggplot2::geom_point(shape=1) + 
                        ggplot2::geom_smooth(method="auto") + ggplot2::theme_minimal()       
                
        }
        
        return(res.plot)
}

