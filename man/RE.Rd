\name{RE}
\alias{RE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the relative expression profiles of any given gene expression set.}
\description{
This function computes the relative expression profiles of any given gene expression set. The relative expression profile is being computed as follows:

\deqn{f_s = ( e_s - e_min ) / ( e_max - e_min )}

where \eqn{e_min} and \eqn{e_max} is the minimum/maximum mean expression level over the  developmental stages s. This linear transformation corresponds to a shift by \eqn{e_min} and a subsequent shrinkage by \eqn{e_max - e_min}. As a result, the relative expression level \eqn{f_s} of developmental stage s with minimum \eqn{e_s} is 0, the relative expression level \eqn{f_s} of the developmental stage s with maximum \eqn{e_s} is 1, and the relative expression levels \eqn{f_s} of all other stages s range between 0 and 1, accordingly.

}
\usage{
RE(ExpressionMatrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ExpressionMatrix}{
a numeric matrix representing a gene expression matrix for which the relative expression profile shall be computed.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a vector containing the relative expression profile of the correspnding data matrix.
}
\references{
Domazet-Loso T and Tautz D. 2010. "A phylogenetically based transcriptome age index mirrors ontogenetic divergence patterns". Nature (468): 815-818.

Quint M et al. 2012. "A transcriptomic hourglass in plant embryogenesis". Nature (490): 98-101.

}
\author{
Hajk-Georg Drost}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{REMatrix}}, \code{\link{PlotRE}}
}
\examples{

# read standard phylotranscriptomics data
data(PhyloExpressionSetExample)

# relative expression profile of PS1 genes
RE(PhyloExpressionSetExample[ which(PhyloExpressionSetExample[ , 1] == 1), 3:9 ])


}

