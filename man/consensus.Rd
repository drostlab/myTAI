% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gatai_convergence_plots.R
\name{consensus}
\alias{consensus}
\title{Calculate Consensus Gene Set}
\usage{
consensus(x, p = 0.5)
}
\arguments{
\item{x}{List of gene sets from different GATAI runs}

\item{p}{Frequency threshold (default: 0.5)}
}
\value{
Character vector of consensus genes
}
\description{
Calculate consensus genes from multiple GATAI runs based on
frequency threshold.
}
\details{
This function identifies genes that appear in at least p proportion of
the input gene sets, providing a consensus set of genes across multiple
GATAI runs.
}
\examples{
# Calculate consensus from multiple runs
# consensus_genes <- consensus(gatai_runs, p = 0.5)

}
\keyword{internal}
