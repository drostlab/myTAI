% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{rcpp_tei_parallel}
\alias{rcpp_tei_parallel}
\title{rcpp_tei_parallel}
\usage{
rcpp_tei_parallel(expression, ps, ncores = 1L)
}
\arguments{
\item{expression}{ExpressionSet as sparseMatrix}

\item{ps}{named Phylostratum}

\item{ncores}{number of cores}
}
\value{
list
}
\description{
computes the phylogenetically based
transcriptome evolutionary index (TEI)
}
\examples{
## load example sequence data
data("PhyloExpressionSetExample", package="myTAI")
spmat <- as(data.matrix(PhyloExpressionSetExample[,-c(1,2)]), "sparseMatrix")
rownames(spmat) <- PhyloExpressionSetExample$GeneID
ps <- setNames(PhyloExpressionSetExample$Phylostratum, PhyloExpressionSetExample$GeneID)
rcpp_tei_parallel(spmat, ps)
}
\author{
Kristian K Ullrich
}
