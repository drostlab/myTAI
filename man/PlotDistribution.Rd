\name{PlotDistribution}
\alias{PlotDistribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to plot the frequency distribution of genes within the corresponding \emph{phylostratigraphic map} or \emph{divergence map}.}
\description{
This function plots the frequency distribution of genes within the corresponding \emph{phylostratigraphic map} or \emph{divergence map} and can be used to fastly visualize the PS or DS distribution of a given phylostratum vector or divergence-stratum vector.

}
\usage{
PlotDistribution(PhyloExpressionSet, plotText = TRUE, as.ratio = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{PhyloExpressionSet}{
a standard PhyloExpressionSet object.
}
  \item{plotText}{
a boolean value specifying whether the total number of genes belonging to a corresponding phylostratum or divergence-stratum class
shall be plotted above each bar of the barplot.}

\item{as.ratio}{
a boolean value specifying whether the relative frequencies instead of absolute frequencies shall be plotted.  
}
  \item{\dots}{
default plot parameters.
}
}
\details{
The frequency distribution of all genes or a subset of genes might be of interest for
subsequent analyses.

For Example:

Filtering genes using gene cluster algorithms can result in different groups (classes) of genes.
For each gene group the phylostratum or divergence-stratum distribution can be visualized using this function
and can be compared between different groups.

This analysis allows to compare different gene expression profiles (or gene groups in general) based
on their evolutionary origins or evolutionary relationships.

}
\value{
a barplot showing the phylostratum distribution or divergence-stratum distribution of a given numeric vector containing PS or DS values.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Hajk-Georg Drost
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{PlotPhyloMap}}
}
\examples{
# load PhyloExpressionSet
data(PhyloExpressionSetExample)

# plot the phylostratum distribution of a PhyloExpressionSet
PlotDistribution(PhyloExpressionSetExample, plotText = TRUE)

# plot the relative frequency distribution of a PhyloExpressionSet
PlotDistribution(PhyloExpressionSetExample, plotText = TRUE, as.ratio = TRUE)


# a example for visualizing the PS distribution for a subset of genes
# plot the relative frequency distribution of a PhyloExpressionSet
PlotDistribution(subset(PhyloExpressionSetExample,Phylostratum == 5), plotText = TRUE, as.ratio = TRUE)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
