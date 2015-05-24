#' @title An Example DivergenceExpressionSet Data Set
#' @docType data
#' @description A standard DivergenceExpressionSet is a \code{\link{data.frame}} consisting of a standardized sequence of columns
#' to store the age information for each gene and its corresponding gene expression profile.
#'
#' The standard is defined as follows:
#'        
#' Divergencestratum | GeneID | Expression-level 1 | ... | Expression-level N
#' @details 
#' 
#' This example DivergenceExpressionSet dataset covers 7 developmental stages of Arabidopsis thaliana embryo development.
#' The initial gene expression dataset was published by Xiang et al., 2011 (see references section) 
#' and was then used by Quint et al., 2012 (see references section) to assign sequence divergence values (divergence strata) to each gene expression profile.
#' @return a standard DivergenceExpressionSet object.
#' @references 
#' 
#' Quint M et al. 2012. "A transcriptomic hourglass in plant embryogenesis". Nature (490): 98-101.
#' Supplementary Table 2 : http://www.nature.com/nature/journal/v490/n7418/full/nature11394.html
#' 
#' Xiang D et al. 2011. "Genome-Wide Analysis Reveals Gene Expression and Metabolic Network Dynamics during Embryo Development in Arabidopsis". Plant Physiology (156): 346-356.
#' Supplemental Table 1 : http://www.plantphysiol.org/content/156/1/346/suppl/DC1
#' 
#' 
#' @author Hajk-Georg Drost
#' @seealso \code{\link{PhyloExpressionSetExample}}
#' @name DivergenceExpressionSetExample
#' @source http://www.plantphysiol.org/content/156/1/346/suppl/DC1
NULL



#' @title An Example PhyloExpressionSet Data Set
#' @docType data
#' @description A standard PhyloExpressionSet is a \code{\link{data.frame}} consisting of a standardized sequence of columns
#' to store the age information for each gene and its corresponding gene expression profile.
#'
#' The standard is defined as follows:
#'        
#' Phylostratum | GeneID | Expression-level 1 | ... | Expression-level N
#' @details This dataset covers 7 developmental stages of Arabidopsis thaliana embryo development.
#' The initial gene expression dataset was published by Xiang et al., 2011 (see references section) and was then
#' used by Quint et al., 2012 (see references section) to assign evolutionary ages to each gene expression profile.
#' @return a standard PhyloExpressionSet object.
#' @references 
#' 
#' Quint M et al. 2012. "A transcriptomic hourglass in plant embryogenesis". Nature (490): 98-101.
#' Supplementary Table 2 : http://www.nature.com/nature/journal/v490/n7418/full/nature11394.html
#' 
#' Xiang D et al. 2011. "Genome-Wide Analysis Reveals Gene Expression and Metabolic Network Dynamics during Embryo Development in Arabidopsis". Plant Physiology (156): 346-356.
#' Supplemental Table 1 : http://www.plantphysiol.org/content/156/1/346/suppl/DC1
#' @author Hajk-Georg Drost
#' @seealso \code{\link{DivergenceExpressionSetExample}}
#' @name PhyloExpressionSetExample
#' @source http://www.plantphysiol.org/content/156/1/346/suppl/DC1
NULL



