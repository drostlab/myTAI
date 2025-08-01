% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_relative_expression.R
\name{plot_relative_expression_bar}
\alias{plot_relative_expression_bar}
\title{Plot Mean Relative Expression Levels as Barplot}
\usage{
plot_relative_expression_bar(phyex_set, groups, p_adjust_method = NULL, ...)
}
\arguments{
\item{phyex_set}{A PhyloExpressionSet object (BulkPhyloExpressionSet or ScPhyloExpressionSet).}

\item{groups}{A list of integer vectors specifying age categories (e.g., phylostrata) for each group (2+ groups).}

\item{p_adjust_method}{P-value adjustment for multiple testing.}

\item{...}{Further arguments passed to ggplot2 geoms.}
}
\value{
ggplot2 object.
}
\description{
Plots mean relative expression levels for age category groups using a PhyloExpressionSet S7 object, with statistical testing.
}
