\name{TDI}
\alias{TDI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the Transcriptome Divergence Index (TDI).}
\description{
This function computes the sequence distance based transcriptome divergence index (TDI) originally published by
Quint et al., 2012.}
\usage{
TDI(DivergenceExpressionSet)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{DivergenceExpressionSet}{
a standard PhyloExpressionSet or DivergenceExpressionSet object.
}
}
\details{
The TDI measure represents the weighted arithmetic mean (expression levels as
weights for the divergence-stratum value) over all gene divergence categories denoted as \emph{divergence-strata}.



\deqn{TDI_s = \sum (e_is * ds_i) / \sum e_is}



where TDI_s denotes the TDI value in developmental stage s, e_is denotes the gene expression level of gene i in stage s, and ds_i denotes the corresponding divergence-stratum of gene i, \eqn{i = 1,...,N} and N = total number of genes.




Internally the function calls the C++ function \code{\link{cpp_TAI}} to speed up TDI computations.
}
\value{
a numeric vector containing the TDI values for all given developmental stages.
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
\code{\link{TAI}}, \code{\link{cpp_TAI}}
}
\examples{
# reading a standard DivergenceExpressionSet
data(DivergenceExpressionSetExample)

# computing the TDI profile of a given DivergenceExpressionSet object
TDIs <- TDI(DivergenceExpressionSetExample)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
