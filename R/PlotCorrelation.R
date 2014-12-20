#' @title Plot the correlation between Phylostrata and Divergence-Strata
#'  of a given PhyloExpressionSet and DivergenceExpressionSet
#'  @description
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
#' @param main.text a character string specifying the text that shall be written as main. 
#' Default is \code{main.text} = \code{NULL}.
#' @param \dots default plot parameters.
#' @return a jitter-correlation-plot of PS and DS correlation.
#' @references  Quint M et al. 2012. "A transcriptomic hourglass in plant embryogenesis". Nature (490): 98-101.
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
#'                 method = "pearson", linearModel = TRUE, 
#'                 main.text = "Pearson's ")
#' 
#' 
#' 
#' @export

PlotCorrelation <- function(PhyloExpressionSet,DivergenceExpressionSet,
                            method = "pearson",linearModel = FALSE, 
                            main.text = NULL,...)
{
        
        is.ExpressionSet(PhyloExpressionSet)
        is.ExpressionSet(DivergenceExpressionSet)
        
        
        if(!is.element(method, c("pearson", "kendall", "spearman")))
                stop("Please choose a correlation method that is supported by this function.")
        
        colnames(PhyloExpressionSet)[2] <- "GeneID"
        colnames(DivergenceExpressionSet)[2] <- "GeneID"
        
        # convert ids to lower case
        PhyloExpressionSet[ , "GeneID"] <- tolower(PhyloExpressionSet[ , "GeneID"])
        DivergenceExpressionSet[ , "GeneID"] <- tolower(DivergenceExpressionSet[ , "GeneID"])
        
        PS_DS.Subset <- merge(PhyloExpressionSet[ , 1:2], DivergenceExpressionSet[ , 1:2],by = "GeneID")
        
        CorrelationCoefficient <- cor(PS_DS.Subset[ , 2],PS_DS.Subset[ , 3],method = method)
        CorrCoeffasCharacter <- as.character(round(CorrelationCoefficient,3))
        
        nrows <- dim(PS_DS.Subset)[1]
        PS <- vector(mode = "numeric", length = nrows)
        DS <- vector(mode = "numeric", length = nrows)
        
        PS <- jitter(PS_DS.Subset[ , 2],1.5)
        DS <- jitter(PS_DS.Subset[ , 3],1.5)
        
        if(!is.null(main.text))
                plot(PS,DS,main = paste0(main.text,method," = ",CorrCoeffasCharacter), ...)
        
        if(is.null(main.text))
                plot(PS,DS,main = paste0(method," = ",CorrCoeffasCharacter), ...)
        
        if(linearModel == TRUE)
                abline(lm(PS_DS.Subset[ , 3]~PS_DS.Subset[ , 2]),lwd = 5,col = "red")
        
}

