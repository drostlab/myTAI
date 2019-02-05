#' @title Plot Transcriptome Index using bootstrapping and confidence intervals
#' @description 
#' Function to plot and compare the confidence intervals of Transcriptome Index between transformed and non-transformed expression data 
#' by using bootstrapping appraoches instead of permutation tests used in \code{\link{PlotSignature}}.
#' @details 
#' This function can be used to check potential outliers (e.g. a few exramly highly expressed genes) in transcriptome. 
#' Since Transcriptome Index is weigthed mean, it could be easily affectd by outliers. So, we have to check potential outliers in the transcriptome data.
#' Because log or sqrt trandformation can alleviate the effect of outliers, if there are some outliers, we could see the confidence intervals (genetated by bootstrapping) 
#' from non-trandformed expression data are much higher and more variable than from log or sqrt trandformed expression data.
#' In order to compare the range of confidence intervals in the same scale, we plotted the ratio of upper to lower confidence interval boundary across development.
#' @param ExpressionSet a standard PhyloExpressionSet, DivergenceExpressionSet or PolymorphismsExpressionSet object.
#' @param measure type of transcriptome index that shall be computed. E.g. measure = "TAI" (Transcriptome Age Index), measure = "TDI" (Transcriptome Divergence Index), measure = "TPI" (Transcriptome Polymorphism Index).
#' @param nbootstraps number of independent bootstraps.
#' @author Jialin Liu
#' @examples 
#' data("PhyloExpressionSetExample")
#' PlotCIRatio(PhyloExpressionSetExample,"TAI",5)
#' @seealso \code{\link{PlotSignature}}
#' @export

PlotCIRatio <- function(ExpressionSet, measure, nbootstraps) {
        
        if (!is.element(measure, c("TAI", "TDI", "TPI")))
                stop(
                        "Measure '",
                        measure,
                        "' is not available for this function. Please specify a measure supporting by this function.",
                        call. = FALSE
                )
        
        ExpressionSet <- as.data.frame(ExpressionSet)
        names(ExpressionSet)[2] <- "GeneID"
        bootID <-
                replicate(nbootstraps,
                          sample(ExpressionSet$GeneID, replace = TRUE))
        rawIndex <- c()
        logIndex <- c()
        sqrtIndex <- c()
        for (i in 1:nbootstraps) {
                tempID <- data.frame(bootID[ , i])
                names(tempID) <- "GeneID"
                tempData <-
                        merge(ExpressionSet, tempID, by = "GeneID")
                tempData <- tempData[c(2, 1, 3:ncol(ExpressionSet))]
                if (measure == "TAI") {
                        ## raw expression value
                        rawIndex <- rbind(rawIndex, TAI(tempData))
                        ## log2 expression value
                        tempData_log2 <- tf(tempData, function(x) log2(x + 1))
                        logIndex <-
                                rbind(logIndex, TAI(tempData_log2))
                        ## square root expression value
                        tempData_sqrt <- tf(tempData, sqrt)
                        sqrtIndex <-
                                rbind(sqrtIndex, TAI(tempData_sqrt))
                }
                if (measure == "TDI") {
                        ## raw expression value
                        rawIndex <- rbind(rawIndex, TDI(tempData))
                        ## log2 expression value
                        tempData_log2 <- tf(tempData, function(x) log2(x + 1))
                        logIndex <-
                                rbind(logIndex, TDI(tempData_log2))
                        ## square root expression value
                        tempData_sqrt <- tf(tempData, sqrt)
                        sqrtIndex <-
                                rbind(sqrtIndex, TDI(tempData_sqrt))
                }
        }
        ## interval ratio
        rawRatio <-
                apply(rawIndex, 2, function(x)
                        stats::quantile(x,  probs = 0.975)) /
                apply(rawIndex, 2, function(x)
                    stats::quantile(x,  probs = 0.025))
        logRatio <-
                apply(logIndex, 2, function(x)
                    stats::quantile(x,  probs = 0.975)) /
                apply(logIndex, 2, function(x)
                    stats::quantile(x,  probs = 0.025))
        sqrtRatio <-
                apply(sqrtIndex, 2, function(x)
                    stats::quantile(x,  probs = 0.975)) /
                apply(sqrtIndex, 2, function(x)
                    stats::quantile(x,  probs = 0.025))
        ## plot
        graphics::par(mfrow = c(1, 1))
        graphics::par(mar = c(7, 5, 2, 2))
        allRatio <- c(rawRatio, logRatio, sqrtRatio)
        graphics::plot(
                logRatio,
                type = "l",
                col = "black",
                lwd = 5,
                ylim = c(min(allRatio), max(allRatio) + max(allRatio) * 0.05),
                xlab = "",
                ylab = " CI boundary ratio",
                cex.lab = 1.4,
                cex.axis = 1.2,
                xaxt = 'n'
        )
        graphics::lines(sqrtRatio,
              col = "black",
              lty = 2,
              lwd = 5)
        graphics::lines(rawRatio,
              col = "black",
              lty = 3,
              lwd = 5)
        graphics::legend(
                "topleft",
                c(
                        "Without transformation",
                        "Srqt transformation",
                        "Log transformation"
                ),
                lty = c(3, 2, 1),
                col = "black",
                lwd = 2,
                cex = 1.4,
                bty = "n"
        )
        devNames <- colnames(ExpressionSet)[c(-1,-2)]
        for (j in 1:length(devNames)) {
                graphics::axis(
                        side = 1,
                        at = j,
                        labels = devNames[j],
                        las = 2,
                        cex.axis = 1.2
                ) # Add development stages as labels
        }
}

