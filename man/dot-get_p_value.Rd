% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stat_test_result.R
\name{.get_p_value}
\alias{.get_p_value}
\title{Calculate P-Value from Distribution}
\usage{
.get_p_value(
  cdf,
  test_stat,
  params,
  alternative = c("two-sided", "less", "greater")
)
}
\arguments{
\item{cdf}{Cumulative distribution function}

\item{test_stat}{Numeric test statistic value}

\item{params}{List of distribution parameters}

\item{alternative}{Character string specifying alternative hypothesis 
("two-sided", "less", "greater")}
}
\value{
Numeric p-value
}
\description{
Internal function to calculate p-values from cumulative distribution functions
based on test statistics and alternative hypothesis specifications.
}
\details{
This function calculates p-values using the appropriate tail(s) of the distribution:
- "greater": Uses upper tail (1 - CDF)
- "less": Uses lower tail (CDF)  
- "two-sided": Uses 2 * minimum of both tails
}
\examples{
# Calculate p-value (internal use)
# pval <- .get_p_value(pnorm, 1.96, list(mean=0, sd=1), "two-sided")

}
\keyword{internal}
