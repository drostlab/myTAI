\name{MatchMap}
\alias{MatchMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to match a phylostratigraphic map  or divergence map with an ExpressionSet to obtain a standard PhyloExpressionSet or DivergenceExpressionSet object.}
\description{
This function matches a \emph{phylostratigraphic map} or \emph{divergence map} only storing unique gene ids with an ExpressionSet also storing only unique gene ids.

}
\usage{
MatchMap(Map, ExpressionMatrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Map}{
a standard \emph{phylostratigraphic map} or \emph{divergence map} object.
}
  \item{ExpressionMatrix}{
a standard ExpressionMatrix object.
}
}
\details{
In phylotranscriptomics analyses two major techniques are performed to
obtain standard \emph{phylostratigraphic map} or \emph{divergence map} objects.

To obtain a \emph{phylostratigraphic map}, \emph{phylostratigraphy} (Domazet-Lošo et al., 2007) has to be performed. To obtain a \emph{sequence divergence map}, orthologous gene detection and Ka/Ks computations (Quint et al., 2012; Drost et al., 2014) have to be performed.

The resulting standard \emph{phylostratigraphic map} or \emph{divergence map} objects consist of 2 colums storing the phylostratum assignment or divergence stratum assignment of a given gene in column one, and the corresponding gene id of that gene on columns two.

A standard ExpressionMatrix is a common gene expression matrix storing the gene ids or probe ids in the first column, and all experiments/stages/replicates in the following columns.

The \emph{MatchMap} function takes both standard datasets: \emph{Map} and \emph{ExpressionMatrix} to merge both data sets to obtain a standard PhyloExpressionSet or DivergenceExpressionSet object.

This procedure is analogous to \code{\link{merge}}, but is customized to the \emph{phylostratigraphic map}, \emph{divergence map}, and \emph{ExpressionMatrix} standards to allow a faster and more intuitive usage. 
}
\value{
a standard PhyloExpressionSet or DivergenceExpressionSet object.
}
\references{
Domazet-Loso T, Brajković J, Tautz D (2007) A phylostratigraphy approach to uncover the genomic history of major adaptations in metazoan lineages. Trends Genet. 23: 533-9.

Domazet-Loso T, Tautz D (2010) A phylogenetically based transcriptome age index mirrors ontogenetic divergence patterns. Nature 468: 815-8.

Quint M., Drost H.G., Gabel A., Ullrich K.K., Boenn M., Grosse I. (2012) A transcriptomic hourglass in plant embryogenesis. Nature 490: 98-101.

Drost et al. 2014, Active maintenance of phylotranscriptomic hourglass patterns in animal and plant embryogenesis.
}
\author{
Hajk-Georg Drost
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

# load a standard PhyloExpressionSet
data(PhyloExpressionSetExample)

# in a standard PhyloExpressionSet, column one and column two denote a standard phylostratigraphic map
PhyloMap <- PhyloExpressionSetExample[ , 1:2]

# look at the phylostratigraphic map standard
head(PhyloMap)

# in a standard PhyloExpressionSet, column two combined with column 3 - N denote a standard ExpressionMatrix
ExpressionMatrixExample <- PhyloExpressionSetExample[ , c(2,3:9)]

# these two data sets shall illustrate an example phylostratigraphic map that is returned
# by a standard phylostratigraphy run, and a expression set that is the result of expression data analysis (background correction, normalization, ...)

# now we can use the MatchMap function to merge both data sets
# to obtain a standard PhyloExpressionSet

PES <- MatchMap(PhyloMap, ExpressionMatrixExample)

# note that the function returns a head() of the matched gene ids to enable
# the user to find potential mis-matches

# the entire procedure is analogous to merge() with two data sets sharing the same gene ids as column (primary key)
PES_merge <- merge(PhyloMap, ExpressionMatrixExample)



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
