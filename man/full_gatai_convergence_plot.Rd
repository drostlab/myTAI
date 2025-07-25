% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gatai_convergence_plots.R
\name{full_gatai_convergence_plot}
\alias{full_gatai_convergence_plot}
\title{Create Full GATAI Convergence Plot}
\usage{
full_gatai_convergence_plot(phyex_set, runs, p = 0.5, ps = c(0.25, 0.5, 0.75))
}
\arguments{
\item{phyex_set}{A PhyloExpressionSet object}

\item{runs}{List of GATAI run results}

\item{p}{Consensus threshold for petal plot (default: 0.5)}

\item{ps}{Vector of consensus thresholds for convergence plots (default: c(0.25, 0.5, 0.75))}
}
\value{
A cowplot grid showing convergence analysis
}
\description{
Create a comprehensive plot showing GATAI convergence across multiple
runs and thresholds.
}
\details{
This function creates a comprehensive visualization of GATAI convergence
including consensus set sizes, p-values, threshold comparisons, and
gene removal patterns across multiple runs.
}
\examples{
# Create full convergence plot
# conv_plot <- full_gatai_convergence_plot(phyex_set, gatai_runs)

}
\keyword{internal}
