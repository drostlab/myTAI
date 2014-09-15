\name{PlotCorrelation}
\alias{PlotCorrelation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to plot the correlation between phylostratum values and divergence-stratum values
of a given PhyloExpressionSet and DivergenceExpressionSet.}
\description{
This function plots the correlation coefficient between phylostratum values and divergence-stratum values of a given PhyloExpressionSet and DivergenceExpressionSet.

This function can be used to test whether a given PS distribution and DS distribution are linear correlated so that the independence of PS and DS can be assumed for subsequent analyses (Quint et al., 2012).

}
\usage{
PlotCorrelation(PhyloExpressionSet, 
                    DivergenceExpressionSet, method = "pearson",
                    linearModel=F, main.text="", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PhyloExpressionSet}{
a standard PhyloExpressionSet object.
}
  \item{DivergenceExpressionSet}{
a standard DivergenceExpressionSet object.
}
  \item{method}{
a character string specifying the correlation method to cbe used, e.g. 
"pearson", "kendall", "spearman".  
 }

  \item{linearModel}{
a boolean value specifying whether a linear model should be fitted to the data and furthermore, should be visualized in the corresponding plot.
}


\item{main.text}{

a character string specifying the text that shall be written as main.

}
  \item{\dots}{
default plot parameters.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a jitter-plot of PS and DS correlation.
}
\references{
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
\code{\link{cor}}
}
\examples{

# read standard phylotranscriptomics data
data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

# plot the PS and DS correlation
PlotCorrelation(PhyloExpressionSetExample, DivergenceExpressionSetExample, method = "kendall", linearModel = TRUE, main.text = "Kendall's ")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
