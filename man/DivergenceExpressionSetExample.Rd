\name{DivergenceExpressionSetExample}
\docType{data}
\alias{DivergenceExpressionSetExample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A example data set for a standard DivergenceExpressionSet.
}
\description{
A standard DivergenceExpressionSet is a \code{\link{data.frame}} consisting of a standardized sequence of columns
to store the age information for each gene and its corresponding gene expression profile.

The standard is defined as follows:

Divergencestratum | GeneID | Expression-level 1 | ... | Expression-level N

}
\usage{
DivergenceExpressionSetExample
}
%- maybe also 'usage' for other objects documented here.
\arguments{
}
\details{
This dataset covers 7 developmental stages of Arabidopsis thaliana embryo development.
The initial gene expression dataset was published by Xiang et al., 2011 and was then
used by Quint et al., 2012 to assign evolutionary ages to each gene expression profile.
}
\value{
a standard DivergenceExpressionSet object.
}
\references{
Quint M et al. 2012. "A transcriptomic hourglass in plant embryogenesis". Nature (490): 98-101.

Xiang D et al. 2011. "Genome-Wide Analysis Reveals Gene Expression and Metabolic Network Dynamics during Embryo Development in Arabidopsis". Plant Physiology (156): 346-356.
}
\author{
Hajk-Georg Drost
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{PhyloExpressionSetExample}}
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
