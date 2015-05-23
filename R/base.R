
pow <- function(x,power)
{
  return(x^power)
}

CollapseFromTo <- function(x,from,to,FUN){  f <- match.fun(FUN); return(apply(x[ , from:to], 1 , f)) }

#' @title Match a Phylostratigraphic Map or Divergence Map with a ExpressionMatrix
#' @details This function matches a \emph{Phylostratigraphic Map} or \emph{Divergence Map} only storing unique gene ids with a ExpressionMatrix
#' also storing only unique gene ids.
#' @param Map a standard \emph{Phylostratigraphic Map} or \emph{Divergence Map} object.
#' @param ExpressionMatrix  a standard ExpressionMatrix object.
#' @param remove.duplicates a logical value indicating whether duplicate gene ids should be removed from the data set.
#' @param accumulate an accumulation function such as \code{mean()}, \code{median()}, or \code{min()}
#' to accumulate multiple expression levels that map to the same unique gene id present in the \code{ExpressionMatrix}.
#' @details
#' 
#' In phylotranscriptomics analyses two major techniques are performed to
#' obtain standard \emph{Phylostratigraphic map} or \emph{Divergence map} objects.
#' 
#' To obtain a \emph{Phylostratigraphic Map}, \emph{Phylostratigraphy} (Domazet-Loso et al., 2007) has to be performed. To obtain a \emph{Divergence Map}, 
#' orthologous gene detection, Ka/Ks computations, and decilation (Quint et al., 2012; Drost et al., 2015) have to be performed.
#' 
#' The resulting standard \emph{Phylostratigraphic Map} or \emph{Divergence Map} objects consist of 2 colums storing the phylostratum assignment 
#' or divergence stratum assignment of a given gene in column one, and the corresponding gene id of that gene on columns two.
#' 
#' A standard ExpressionMatrix is a common gene expression matrix storing the gene ids or probe ids in the first column, and all experiments/stages/replicates in the following columns.
#' 
#' The \emph{MatchMap} function takes both standard datasets: \emph{Map} and \emph{ExpressionMatrix} to merge both data sets to obtain a standard PhyloExpressionSet or DivergenceExpressionSet object.
#' 
#' This procedure is analogous to \code{\link{merge}}, but is customized to the \emph{Phylostratigraphic Map}, \emph{Divergence Map}, and \emph{ExpressionMatrix} standards to allow a faster and more intuitive usage. 
#' 
#' In case you work with an ExpressionMatrix that stores multiple expression levels for a unique gene id, you
#' can specify the \code{accumulation} argument to accumulate these multiple expression levels to obtain
#' one expression level for one unique gene.
#' 
#' 
#' @return a standard PhyloExpressionSet or DivergenceExpressionSet object. 
#' @references 
#' 
#'   Domazet-Loso T, Brajkovic J, Tautz D (2007) A phylostratigraphy approach to uncover the genomic history of major adaptations in metazoan lineages. Trends Genet. 23: 533-9.
#' 
#'   Domazet-Loso T, Tautz D (2010) A phylogenetically based transcriptome age index mirrors ontogenetic divergence patterns. Nature 468: 815-8.
#' 
#'   Quint M., Drost H.G., Gabel A., Ullrich K.K., Boenn M., Grosse I. (2012) A transcriptomic hourglass in plant embryogenesis. Nature 490: 98-101.
#' 
#'   Drost HG et al. (2015) \emph{Evidence for Active Maintenance of Phylotranscriptomic Hourglass Patterns in Animal and Plant Embryogenesis}. Mol Biol Evol. 32 (5): 1221-1231 doi:10.1093/molbev/msv012.
#' 
#' @author Hajk-Georg Drost
#' @examples
#'         
#' # load a standard PhyloExpressionSet
#' data(PhyloExpressionSetExample)
#'         
#' # in a standard PhyloExpressionSet, 
#' # column one and column two denote a standard 
#' # phylostratigraphic map
#' PhyloMap <- PhyloExpressionSetExample[ , 1:2]
#'         
#' # look at the phylostratigraphic map standard
#' head(PhyloMap)
#'         
#' # in a standard PhyloExpressionSet, column two combined 
#' # with column 3 - N denote a standard ExpressionMatrix
#' ExpressionMatrixExample <- PhyloExpressionSetExample[ , c(2,3:9)]
#'         
#' # these two data sets shall illustrate an example 
#' # phylostratigraphic map that is returned
#' # by a standard phylostratigraphy run, and a expression set 
#' # that is the result of expression data analysis 
#' # (background correction, normalization, ...)
#'         
#' # now we can use the MatchMap function to merge both data sets
#' # to obtain a standard PhyloExpressionSet
#'         
#' PES <- MatchMap(PhyloMap, ExpressionMatrixExample)
#'         
#' # note that the function returns a head() 
#' # of the matched gene ids to enable
#' # the user to find potential mis-matches
#'         
#' # the entire procedure is analogous to merge() 
#' # with two data sets sharing the same gene ids 
#' # as column (primary key)
#' PES_merge <- merge(PhyloMap, ExpressionMatrixExample)
#'         
#'         
#'         
#'@export
MatchMap <- function(Map,ExpressionMatrix, remove.duplicates = FALSE, accumulate = NULL)
{
  
  names(ExpressionMatrix)[1] <- "GeneID"
  names(Map)[2] <- "GeneID"
  ExpressionMatrix[ , "GeneID"] <- tolower(ExpressionMatrix[ , "GeneID"])
  Map[ , "GeneID"] <- tolower(Map[ , "GeneID"])
  GeneID <- NULL
  
  if(remove.duplicates)
          Map <- Map[-which(duplicated(Map[ , "GeneID"])) , ]
  
  if(any(duplicated(Map[ , "GeneID"])))
          stop("You have duplicate Gene IDs in your Map. Please enter only unique Gene IDs.")
  
  if(!is.null(accumulate)){
        
          acc_fun <- match.fun(accumulate)
          ExpressionMatrix <- dplyr::summarise_each(dplyr::group_by(ExpressionMatrix, GeneID), dplyr::funs(acc_fun))
           
  }
          
  if(any(duplicated(ExpressionMatrix[ , "GeneID"])))
          stop("You have duplicate Gene IDs in your ExpressionMatrix. Please enter only unique Gene IDs, or specify the 'accumulate' argument.")
  
  joined_ExpressionMatrix <- dplyr::semi_join(ExpressionMatrix, Map, by = "GeneID")
  
  res_tbl <- merge(joined_ExpressionMatrix,Map, by = "GeneID")
  
  if(!any(duplicated(res_tbl[ , "GeneID"]))){
          
          return(res_tbl[ , c(ncol(res_tbl),1:(ncol(res_tbl)-1))])
          
  } else {
          stop("Something went wrong with matching Map and ExpressionMatrix! Please check for duplicate entries!")
  }
  
}


#' @title Compute TAI or TDI Profiles Omitting a Given Gene
#' @description For each gene i, exclude the corresponding gene i from the global
#'  PhyloExpressionSet or DivergenceExpressionSet and compute the \code{\link{TAI}} or \code{\link{TDI}} 
#'  profile for the corresponding global PhyloExpressionSet or DivergenceExpressionSet
#'  with excluded gene i. 
#'  
#'  This procedure results in a TAI or TDI profile Matrix storing the TAI or TDI profile for each omitted gene i.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet object.        
#' @return a numeric matrix storing TAI or TDI profile for each omitted gene i.
#' @author Hajk-Georg Drost
#' @examples
#' 
#' # read standard phylotranscriptomics data
#' data(PhyloExpressionSetExample)
#' data(DivergenceExpressionSetExample)
#'
#' # example PhyloExpressionSet
#' omMatrix_ps <- omitMatrix(PhyloExpressionSetExample)
#'
#' # example DivergenceExpressionSet
#' omMatrix_ds <- omitMatrix(DivergenceExpressionSetExample)
#' 
#' 
#' @export
omitMatrix <- function(ExpressionSet)
{
  
  is.ExpressionSet(ExpressionSet)
  
  ncols <- dim(ExpressionSet)[2]
  
  oMatrix <- matrix(NA_real_, ncol = (ncols - 2), nrow = nrow(ExpressionSet))
  oMatrix <- cpp_omitMatrix(as.matrix(ExpressionSet[ , 3:ncols]),as.vector(ExpressionSet[ , 1]))
  
  colnames(oMatrix) <- names(ExpressionSet)[3:ncols]
  rownames(oMatrix) <- paste0("(-) ",ExpressionSet[ , 2])
          
  return(oMatrix)
  
}


quant <- function(ExpressionMatrix,quantile = 0.9)
{
        nCols <- dim(ExpressionMatrix)[2]
        threshold <- vector(mode = "numeric",length = nCols)
        
        for(i in 1:nCols){
                threshold[i] <- as.numeric(quantile(ExpressionMatrix[ , i],probs = quantile))
        }
        return (threshold)
}

std.error <- function(x)
{
        if(is.numeric(x)){
                return(cpp_std_error(as.vector(x)))
        }
        else{
                stop("Please enter a numeric vector.")
        }
}

geom.mean <- function(x)
{
        
        if(is.numeric(x)){
                return(cpp_geom_mean(as.vector(x)))
        }
        
        else{
                stop("Please enter a numeric vector.")
        }
        
}

harmonic.mean <- function(x)
{
        
        N <- length(x)
        ratio <- 1/x
        harm <- 1 / ((1/N) * sum(ratio))
        return(harm)
        
}


re.colors <- function(n)
{
        
        colos <- c("black","red","green","brown","darkmagenta",
                   "blue","darkred","darkblue","darkgreen", "orange",
                   "azure4","gold4","greenyellow","hotpink4",
                   "mediumorchid3","mediumorchid3","peachpuff4",
                   "hotpink","lightgoldenrod", "peru", "slateblue3", "yellow4", "yellowgreen")
        
        return(colos[1:n])
        
}


#' @title Color palette for barplots
#' @description A nice color palette for barplots with several bars.
#' @param n the number of colors to be in the palette. 
#' @return a character vector containing different color names that can be used for barplots.
#' @details This function can be used to select colors for bar plots. 
#' @author Hajk-Georg Drost
#' @examples
#' 
#' # get 5 different colors for 5 different bars
#' barplot_colors <- bar.colors(5)
#' @export
bar.colors <- function(n)
{
        
        colos <- c("black","gray76","gray43","navy","lightskyblue3","palegreen4","seagreen1","lavenderblush1")
        return(colos[1:n])
        
}


#' @title Test ExpressionSet Standard
#' @description This function tests whether a given ExpressionSet follows the pre-defined PhyloExpressionSet or DivergenceExpressionSet standard.
#' @param ExpressionSet a standard PhyloExpressionSet or DivergenceExpressionSet that shall be tested for format validity.
#' @author Hajk-Georg Drost
#' @examples
#' 
#' # read example PhyloExpressionSet
#' data(PhyloExpressionSetExample)
#' 
#' is.ExpressionSet(PhyloExpressionSetExample)
#' 
#' @export
is.ExpressionSet <- function(ExpressionSet){
        
        ncols <- dim(ExpressionSet)[2]
        
        d.f_bool <- is.data.frame(ExpressionSet)
        age.vector_bool <- is.numeric(ExpressionSet[ , 1])
        gene.vector_bool <- ifelse(is.factor(ExpressionSet[ , 2]),is.character(levels(ExpressionSet[ , 2])),is.character(ExpressionSet[ , 2]))
        expression.matrix_bool <- all(sapply(ExpressionSet[ , 3:ncols], is.numeric))
        any.NA.values_bool <- !any(is.na(ExpressionSet))
        
        
        if(all(c(d.f_bool,age.vector_bool,gene.vector_bool,expression.matrix_bool,any.NA.values_bool))){
                return(TRUE)
        }
        
        else{
                stop("The present input object does not fulfill the ExpressionSet standard.")
        }
}




