# Transcriptome age index
TAI <- function(PhyloExpressionSet)
{
  
  is.ExpressionSet(PhyloExpressionSet)
  
  nCols <- dim(PhyloExpressionSet)[2]
  ExpressionMatrix <- PhyloExpressionSet[ , 3:nCols]
  Phylostratum <- PhyloExpressionSet[ , 1]
  TAIProfile <- vector(mode = "numeric",length = nCols-2)
  
  TAIProfile <- cpp_TAI(as.matrix(ExpressionMatrix),as.vector(Phylostratum))
  names(TAIProfile) <- names(ExpressionMatrix)
  
  return(TAIProfile)
  
}

# Transcriptome divergence index
TDI <- function(DivergenceExpressionSet)
{
  
  is.ExpressionSet(DivergenceExpressionSet)
  
  nCols <- dim(DivergenceExpressionSet)[2]
  ExpressionMatrix <- DivergenceExpressionSet[ , 3:nCols]
  Divergencestratum <- DivergenceExpressionSet[ , 1]
  TDIProfile <- vector(mode = "numeric",length = nCols-2)
  
  
  TDIProfile <- cpp_TAI(as.matrix(ExpressionMatrix),as.vector(Divergencestratum))
  names(TDIProfile) <- names(ExpressionMatrix)
  
  return(TDIProfile)
  
}

# Computing the Partial TAI Values for each gene and each developmental stage
# The sum over all Partial TAI Values = the total TAI value (TAIs)
pMatrix <- function(ExpressionSet)
{
  
  is.ExpressionSet(ExpressionSet)
  
  nCols <- dim(ExpressionSet)[2]
  nRows <- dim(ExpressionSet)[1]
  pTAIMatrix <- matrix(nrow = nRows,ncol = nCols-2)
  
  pTAIMatrix <- cpp_pMatrix(as.matrix(ExpressionSet[ , 3:nCols]),as.vector(ExpressionSet[ , 1]))
  
  colnames(pTAIMatrix) <- names(ExpressionSet)[3:nCols]
  rownames(pTAIMatrix) <- ExpressionSet[ , 2]
  
  return(pTAIMatrix)
  
}



RE <- function(ExpressionMatrix)
{
  mDimensions <- dim(ExpressionMatrix)
  cMeans <- vector(mode = "numeric",length=mDimensions[2])
  cMeans <- colMeans(ExpressionMatrix)
  
  f_max <- max(cMeans)
  f_min <- min(cMeans)
  RE <- (cMeans - f_min) / (f_max - f_min)
  return(RE)
}


REMatrix <- function(ExpressionSet)
{
  
  is.ExpressionSet(ExpressionSet)
  return(age.apply(ExpressionSet = ExpressionSet, RE))
  
}



# Function to plot or return the mean expression levels over all expression levels
# for each developmental stage
MeanExpression <- function(ExpressionSet,plot=TRUE)
{
  
  is.ExpressionSet(ExpressionSet)
  
  nCols <- dim(ExpressionSet)[2]
  StageNames <- names(ExpressionSet)[3:nCols]
  # computing the mean expression levels for each developmental stage
  MeanExpValues <- colMeans(ExpressionSet[ , 3:nCols])
  
  # computing the standard deviation of the expression levels for each developmental stage
  SDofExpVals <- apply(ExpressionSet[ , 3:nCols],2,sd)
  
  ylim_MIN <- min(MeanExpValues) - (min(MeanExpValues)/10)
  ylim_MAX <- max(MeanExpValues) + (max(MeanExpValues)/7)
  
  # plotting the mean relative expression levels +- standard deviation
  if(plot == TRUE){
    
    plot(MeanExpValues,type = "l",pch = 20,lwd = 8,col = "steelblue",
         ylim = c(ylim_MIN,ylim_MAX),xlab = "Ontogeny",
         ylab = "Mean Expression Level",xaxt = "n",cex = 4)
    
    axis(1,col.axis = "black",at = format(seq(1,length(StageNames),by = 1)),
         labels = StageNames,las = 1,cex.lab = 3,cex.axis = 1.4)
    #abline(h=ylim_MIN,col="black",lty="dotted")
    legend(x = "topleft",legend = c("mean expression"),col = "steelblue",bty = "n",lwd = 6,ncol = 1,cex = 1.5)
    
  }
  
  if(plot == FALSE){
    result_list <- vector("list",length = 2)
    result_list <- list(mean_expression = MeanExpValues,sd = SDofExpVals)
    return(result_list)
  }
}


#### Small function to compute the power of a given variable x
pow <- function(x,power)
{
  return(x^power)
}

#' @title Function to compute TAI or TDI profiles based on N permuted PhyloExpressionSets or DivergenceExpressionSets.
#' @description This function computes the TAI for a row permutated PhyloExpressionSet or DivergenceExpressionSet.
#'
#' One can determine the number of permutations. The function then returns a \code{\link{TAI}} or \code{\link{TDI}} matrix holding
#' the \code{\link{TAI}} or \code{\link{TDI}} profiles of the permutated PhyloExpressionSets or DivergenceExpressionSets. This procedure
#' can be used for Test-Statistics based on the TAI or TDI profiles. 
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param permutations a numeric value specifying the number of permutations to be performed.
#' @details The sampled \code{\link{TAI}} or \code{\link{TDI}} matrix samples the phylostratum or divergence-stratum vector of
#' a given PhyloExpressionSet or DivergenceExpressionSet and computes the corresponding TAI or TDI profiles
#' of the randomly assigned phylostrata or divergence-strata. This sampling is then performed N times, yielding N randomly sampled TAI or TDI profiles.
#' This random TAI or TDI profile matrix can furthermore be used to perform statistical tests (such as the \code{\link{FlatLineTest}} or \code{\link{ReductiveHourglassTest}}) based on the significance of TAI or TDI patterns.
#' @return a numeric matrix representing N randomly permuted TAI or TDI profiles.
#' @references Quint M et al. 2012. "A transcriptomic hourglass in plant embryogenesis". Nature (490): 98-101.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{FlatLineTest}}, \code{\link{ReductiveHourglassTest}}
#' @examples \dontrun{
#' 
#' # read standard phylotranscriptomics data
#' data(PhyloExpressionSetExample)
#' data(DivergenceExpressionSetExample)
#' 
#' # example PhyloExpressionSet using 1000 permutations
#' randomTAI.Matrix <- bootMatrix(PhyloExpressionSetExample, permutations = 1000)
#' 
#' # example DivergenceExpressionSet using 1000 permutations
#' randomTDI.Matrix <- bootMatrix(DivergenceExpressionSetExample, permutations = 1000
#' 
#' 
#' }
#' @export
bootMatrix <- function(ExpressionSet,permutations=1000)
{
  
  is.ExpressionSet(ExpressionSet)
  
  nCols <- dim(ExpressionSet)[2]
  bootstrapMatrix <- matrix(NA_real_, permutations, (nCols - 2))
  ExprSet <- as.matrix(ExpressionSet[ , 3:nCols])
  AgeVector <- as.vector(ExpressionSet[ , 1])
  
  bootstrapMatrix <- cpp_bootMatrix(ExprSet, AgeVector, permutations)
  colnames(bootstrapMatrix) <- names(ExpressionSet)[3:nCols]
  rownames(bootstrapMatrix) <- 1:permutations
  
  return(bootstrapMatrix)
  
  
}



#' @title Function to apply analyses on each phylotranscriptomic age class.
#' @description 
#' This function takes a standard PhyloExpressionSet or DivergenceExpressionSet object
#' as argument and performs any given function to the sub expression matrix of each phylotranscriptomic age class.
#' 
#' This function is verx useful to perform any phylostratum or divergence-stratum specific analysis.
#' 
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' @param FUN a function to be performed on the corresponding expression matrix of each phylostratum or divergence-stratum.
#' @param ... additional arguments of FUN.
#' @param as.list a boolean value specifying whether the output format shall be a matrix or a list object.
#' @details This function uses the \code{\link{split}} function to subset the expression matrix into
#' phylostratum specific sub-matrices. Internally using \code{\link{lapply}}, any function can
#' be performed to the sub-matrices. The return value of this function is a numeric matrix storing
#' the return values by \code{FUN} for each phylostratum and each developmental stage s. 
#' Note that the input \code{FUN} must be an function that can be applied to a matrix (e.g., \code{\link{colMeans}} or \code{\link{RE}}). 
#' In case you use an ananymous function you coud use function(x) apply(x , 2 , var).
#' @return Either a numric matrix storing the return values of the applied function for each age class
#' or a numeric list storing the return values of the applied function for each age class in a list.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{split}}, \code{\link{tapply}}, \code{\link{lapply}}, \code{\link{RE}}, \code{\link{REMatrix}}
#' @examples \dontrun{
#' 
#'  data(PhyloExpressionSetExample)
#' # Example 1
#' # get the relative expression profiles for each phylostratum
#' age.apply(PhyloExpressionSetExample, RE)
#'
#' # this is analogous to 
#' REMatrix(PhyloExpressionSetExample)

#' # Example 2
#' # compute the mean expression profiles for each phylostratum
#' age.apply(PhyloExpressionSetExample, colMeans)
#'
#' # Example 3
#' # compute the variance profiles for each phylostratum
#' age.apply(PhyloExpressionSetExample, function(x) apply(x , 2 , var))
#'
#' # Example 4
#' # compute the range for each phylostratum
#' # Note: in this case, the range() function returns 2 values for each phylostratum
#' # and each developmental stage, hence one should use the argument 'as.list = TRUE'
#' # to make sure that the results are returned properly 
#' age.apply(PhyloExpressionSetExample, function(x) apply(x , 2 , range), as.list = TRUE)
#' 
#' }
#' @export
age.apply <- function(ExpressionSet,FUN, ... ,as.list = FALSE)
{
  f <- match.fun(FUN)
  ncols <- dim(ExpressionSet)[2]
  s <- split(ExpressionSet, ExpressionSet[ , 1])
  
  if(as.list == FALSE){
    res <- t(as.data.frame(lapply(s , function(x) f(as.matrix(x[ , 3:ncols]) , ...))))
    rownames(res) <- levels(as.factor(ExpressionSet[ , 1]))
  }
  
  if(as.list == TRUE){
    res <- lapply(s , function(x) f(as.matrix(x[ , 3:ncols]) , ...))
    names(res) <- levels(as.factor(ExpressionSet[ , 1]))
  }
  
  return(res)
}



tf <- function(ExpressionSet, FUN){
  
  is.ExpressionSet(ExpressionSet)
  
  ncols <- dim(ExpressionSet)[2]
  f <- match.fun(FUN)
  
  return(data.frame(ExpressionSet[ , 1:2] , apply(ExpressionSet[ , 3:ncols] , 2 , f)))
  
  
}


MatchMap <- function(Map,ExpressionMatrix)
{
  
  # looking for intersecting genes in both data sets
  intersectingGeneIDs <- intersect(sort(toupper(Map[ , 2])),sort(toupper(ExpressionMatrix[ , 1])))
  # saving the IDs in each data set of the intersecting genes
  XmapGeneNumbers <- match(intersectingGeneIDs,toupper(Map[ , 2]))
  ExpressionSetGeneNumbers <- match(intersectingGeneIDs,toupper(ExpressionMatrix[ , 1]))
  
  # saving a subset which only contains the intersecting genes
  Xmap_Intersected <- Map[XmapGeneNumbers,]
  ExpressionSet_Intersected <- ExpressionMatrix[ExpressionSetGeneNumbers , ]
  
  # putting together the Xmap and Expression Set
  FinalIntersectedExpressionSet <- data.frame(Xmap_Intersected,ExpressionSet_Intersected)
  # print out a test case which shows if the matching pocess has been performed correctly
  print(head(FinalIntersectedExpressionSet,n = 20))
  print(length(intersectingGeneIDs))
  nColsFinal <- dim(FinalIntersectedExpressionSet)[2]
  # return the matched data set
  return(FinalIntersectedExpressionSet[ , c(1, 3:nColsFinal)])
  
}



omitMatrix <- function(ExpressionSet)
{
  
  is.ExpressionSet(ExpressionSet)
  
  ncols <- dim(ExpressionSet)[2]
  return(cpp_omitMatrix(as.matrix(ExpressionSet[ , 3:ncols]),as.vector(ExpressionSet[ , 1])))
  
}

