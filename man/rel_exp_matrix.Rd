% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_relative_expression.R
\name{rel_exp_matrix}
\alias{rel_exp_matrix}
\title{Compute Relative Expression Matrix for PhyloExpressionSet}
\usage{
rel_exp_matrix(phyex_set)
}
\arguments{
\item{phyex_set}{A PhyloExpressionSet object (BulkPhyloExpressionSet or ScPhyloExpressionSet).}
}
\value{
A matrix with age categories as rows and identities as columns, containing relative expression values.
}
\description{
Computes relative expression profiles for all age categories in a PhyloExpressionSet.
}
