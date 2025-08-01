% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_sample_space.R
\name{plot_sample_space}
\alias{plot_sample_space}
\title{Plot Sample Space Visualization}
\usage{
plot_sample_space(
  phyex_set,
  method = c("PCA", "UMAP"),
  colour_by = c("identity", "TXI"),
  seed = 42,
  ...
)
}
\arguments{
\item{phyex_set}{A PhyloExpressionSet object (BulkPhyloExpressionSet or ScPhyloExpressionSet)}

\item{method}{Character string specifying the dimensionality reduction method: 
"PCA" or "UMAP" (default: "PCA")}

\item{colour_by}{Character string specifying what to colour by: "identity" (default), 
"TXI"}

\item{seed}{Integer seed for reproducible UMAP results (default: 42)}

\item{...}{Additional arguments passed to specific methods}
}
\value{
A ggplot2 object showing the sample space visualisation
}
\description{
Create a dimensional reduction plot to visualize sample relationships
in gene expression space using PCA or UMAP.
}
\details{
This function performs log1p transformation on expression data, removes genes with
zero variance, and applies the specified dimensionality reduction method. Samples
are coloured by their group assignments or TAI values.
}
\examples{
# Create PCA plot coloured by identity
# pca_plot <- plot_sample_space(phyex_set, method = "PCA", colour_by = "identity")

# Create UMAP plot coloured by TXI
# umap_plot <- plot_sample_space(phyex_set, method = "UMAP", colour_by = "TXI")

}
