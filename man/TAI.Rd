\name{TAI}
\alias{TAI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the Transcriptome Age Index (TAI).}
\description{
This function computes the phylogenetically based transcriptome age index (TAI)  originally published by Domazet-Loso & Tautz, 2010.}
\usage{
TAI(PhyloExpressionSet)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PhyloExpressionSet}{
a standard PhyloExpressionSet object.
}
}
\details{
The TAI measure represents the weighted arithmetic mean (expression levels as
weights for the phylostratum value) over all evolutionary age categories denoted as \emph{phylostra}.



\deqn{TAI_s = \sum (e_is * ps_i) / \sum e_is}



where TAI_s denotes the TAI value in developmental stage s, e_is denotes the gene expression level of gene i in stage s, and ps_i denotes the corresponding phylostratum of gene i, \eqn{i = 1,...,N} and N = total number of genes.



Internally the function calls the C++ function \code{\link{cpp_TAI}} to speed up TAI computations.
}
\value{
a numeric vector containing the TAI values for all given developmental stages.
}
\references{
Domazet-Loso T. and Tautz D. 2010. "A phylogenetically based transcriptome age index mirrors ontogenetic divergence patterns". Nature (468): 815-818.

Quint M et al. 2012. "A transcriptomic hourglass in plant embryogenesis". Nature (490): 98-101.
}
\author{
Hajk-Georg Drost
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{cpp_TAI}}
}
\examples{
# reading a standard PhyloExpressionSet
data(PhyloExpressionSetExample)

# computing the TAI profile of a given PhyloExpressionSet object
TAIs <- TAI(PhyloExpressionSetExample)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
