% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotSignatureMultiple.R
\name{PlotSignatureMultiple}
\alias{PlotSignatureMultiple}
\title{Plot evolutionary signatures across transcriptomes, for multiple expression sets}
\usage{
PlotSignatureMultiple(
  ExpressionSets,
  set.labels,
  measure = "TAI",
  TestStatistic = "FlatLineTest",
  modules = NULL,
  permutations = 1000,
  p.value = TRUE,
  shaded.area = FALSE,
  xlab = "Ontogeny",
  ylab = "Transcriptome Index",
  main = "",
  legend.title = "Expression Sets",
  lwd = 4,
  alpha = 0.1,
  y.ticks = 3
)
}
\arguments{
\item{ExpressionSets}{a list of PhyloExpressionSet, DivergenceExpressionSet or PolymorphismsExpressionSet objects.}

\item{set.labels}{a character vector of labels, one for each given expression set}

\item{measure}{type of transcriptome index that shall be computed. E.g.
\itemize{
\item \code{measure = "TAI"} (Transcriptome Age Index)
\item \code{measure = "TDI"} (Transcriptome Divergence Index)
\item \code{measure = "TPI"} (Transcriptome Polymorphism Index)
}}

\item{TestStatistic}{a string defining the type of test statistics to be used to quantify the statistical significance the present phylotranscriptomics pattern.
Possible values can be:
\itemize{
\item \code{TestStatistic} = \code{"FlatLineTest"} : Statistical test for the deviation from a flat line
\item \code{TestStatistic} = \code{"ReductiveHourglassTest"} : Statistical test for the existence of a hourglass shape (high-low-high pattern)
\item \code{TestStatistic} = \code{"ReverseHourglassTest"} : Statistical test for the existence of a reverse hourglass pattern (low-high-low pattern)
\item \code{TestStatistic} = \code{"EarlyConservationTest"} : Statistical test for the existence of a early conservation pattern (low-high-high pattern)
\item \code{TestStatistic} = \code{"LateConservationTest"} : Statistical test for the existence of a late conservation pattern (high-high-low pattern)
}}

\item{modules}{a list storing three elements for the \code{\link{ReductiveHourglassTest}}, \code{\link{EarlyConservationTest}}, \code{\link{LateConservationTest}},
or \code{\link{ReverseHourglassTest}}: early, mid, and late.
Each element expects a numeric vector specifying the developmental stages
or experiments that correspond to each module. For example:
\itemize{
\item \code{modules} = \code{list(early = 1:2, mid = 3:5, late = 6:7)} divides a dataset storing seven developmental stages into 3 modules.
}}

\item{permutations}{a numeric value specifying the number of permutations to be performed for the \code{\link{FlatLineTest}}, \code{\link{EarlyConservationTest}}, \code{\link{LateConservationTest}}, \code{\link{ReductiveHourglassTest}} or \code{\link{ReverseHourglassTest}}.}

\item{p.value}{a boolean value specifying whether the p-value of the test statistic shall be printed within the legend, for each expression set.}

\item{shaded.area}{a boolean value specifying whether a shaded area shall
be drawn for the developmental stages defined to be the presumptive phylotypic period.}

\item{xlab}{label of x-axis.}

\item{ylab}{label of y-axis.}

\item{main}{figure title.}

\item{legend.title}{legend title.}

\item{lwd}{line width.}

\item{alpha}{transparency of the shaded area and error ribbon (between [0,1]). Default is \code{alpha = 0.1}.}

\item{y.ticks}{number of ticks on the y-axis. Default is \code{ticks = 3}.}
}
\value{
a ggplot object visualising the transcriptome index of each given
expression set, together with its standard deviation per stage,
obtained by permuting the gene ages.
The profiles are shown on the same axes, so that they can be readily compared.
Optionally, the p-value of each profile, with respect to the choice of statistic,
is shown.
}
\description{
\emph{PlotSignatureMultiple} is used to compare multiple 
transcriptomic index profiles (e.g. TAI) over a common process of interest
(e.g. a developmental process). 
The main use case is visualising how removing different subsets of genes from
the expression set can perturb the transcriptome index signal.
}
\examples{
data(PhyloExpressionSetExample)

# remove top 1\% expressed genes
genes.top_expression <- TopExpressionGenes(PhyloExpressionSetExample, p=.99)
PhyloExpressionSetExample.top_removed <- subset(PhyloExpressionSetExample, 
                                        !(GeneID \%in\% genes.top_expression))

expression_sets = list(PhyloExpressionSetExample, PhyloExpressionSetExample.top_removed)
set_labels = c("100\%", "99\%")

# Flat line test
PlotSignatureMultiple(ExpressionSets = expression_sets,
                      set.labels = set_labels,
                      TestStatistic="FlatLineTest",
                      main = "A. thaliana embryogenesis",
                      legend.title = "Top Expressed Genes Quantile")
                      
# Reductive hourglass test
PlotSignatureMultiple(ExpressionSets = expression_sets,
                      set.labels = set_labels,
                      TestStatistic="ReductiveHourglassTest",
                      main = "A. thaliana embryogenesis",
                      legend.title = "Top Expressed Genes Quantile",
                      modules=list(early=1:2, mid=3:5, late=6:7),
                      shaded.area=TRUE)                                                        
}
\author{
Stefan Manolache
}
