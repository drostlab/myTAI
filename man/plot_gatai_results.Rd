% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/destroy_pattern.R
\name{plot_gatai_results}
\alias{plot_gatai_results}
\title{Plot Comprehensive GATAI Results}
\usage{
plot_gatai_results(
  phyex_set,
  gatai_result,
  conservation_test = stat_flatline_test,
  runs_threshold = 0.5,
  signature_plot_type = c("separate", "combined")
)
}
\arguments{
\item{phyex_set}{A PhyloExpressionSet object containing the original gene expression data.}

\item{gatai_result}{Result list from \code{destroy_pattern()}, containing GATAI analysis output.}

\item{conservation_test}{Function for conservation test (default: \code{stat_flatline_test}).}

\item{runs_threshold}{Threshold for gene removal consistency across runs (default: 0.5).}

\item{signature_plot_type}{Type of signature plot: "separate" for individual plots, "combined" for overlay (default: both options).}
}
\value{
A named list of ggplot/cowplot objects and results:
  \item{signature_plots}{Signature plots before/after GATAI and top variance removal}
  \item{heatmap_plot}{Heatmap of GATAI-removed genes}
  \item{profiles_plot}{Gene expression profiles of GATAI-removed genes}
  \item{profiles_plot_facet}{Faceted gene profiles by strata}
  \item{gene_space_plot}{Gene space plot of GATAI-removed genes}
  \item{mean_var_plot}{Mean-variance plot highlighting GATAI-removed genes}
  \item{strata_plot}{Phylostrata distribution plot (log obs/exp) for GATAI-removed genes}
  \item{null_dist_plot}{Null distribution plot with test statistics and p-values}
  \item{convergence_plots}{GATAI convergence plots (if available)}
}
\description{
Create a suite of plots summarizing the effects of GATAI gene removal on phylotranscriptomic patterns.
}
\details{
This function provides a comprehensive visualization of the impact of GATAI gene removal,
including transcriptomic signature plots, gene expression profiles, heatmaps, mean-variance relationships,
phylostrata distributions, conservation test comparisons, and convergence diagnostics.
}
\examples{
# Run GATAI analysis
# gatai_result <- destroy_pattern(phyex_set, num_runs = 20)
# Plot results
# plots <- plot_gatai_results(phyex_set, gatai_result)
# Print signature plots
# print(plots$signature_plots)

}
\author{
Filipa Martins Costa, Stefan Manolache
}
