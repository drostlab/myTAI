### Function to plot the correlation between phylostrata and divergence-strata
### methods = "pearson" or "kendall" or "spearman"
PlotCorrelation <- function(PhyloExpressionSet,DivergenceExpressionSet,method="pearson",linearModel=F,main.text="",...)
{
        
        is.ExpressionSet(PhyloExpressionSet)
        is.ExpressionSet(DivergenceExpressionSet)
        
        colnames(PhyloExpressionSet)[2] <- "GeneID"
        colnames(DivergenceExpressionSet)[2] <- "GeneID"
        
        PS_DS.Subset <- merge(PhyloExpressionSet[ , 1:2], DivergenceExpressionSet[ , 1:2],by = "GeneID")
        
        CorrelationCoefficient <- cor(PS_DS.Subset[ , 2],PS_DS.Subset[ , 3],method = method)
        CorrCoeffasCharacter <- as.character(round(CorrelationCoefficient,3))
        
        nrows <- dim(PS_DS.Subset)[1]
        PS <- vector(mode = "numeric", length = nrows)
        DS <- vector(mode = "numeric", length = nrows)
        
        PS <- jitter(PS_DS.Subset[ , 2],1.5)
        DS <- jitter(PS_DS.Subset[ , 3],1.5)
        
        plot(PS,DS,main = paste(main.text,method," = ",CorrCoeffasCharacter,sep = ""),...)
        
        if(linearModel == TRUE)
                abline(lm(PS_DS.Subset[ , 3]~PS_DS.Subset[ , 2]),lwd = 5,col = "red")
        
}

