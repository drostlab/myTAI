\name{cpp_bootMatrix}
\alias{cpp_bootMatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Internal C++ function to compute the bootstrapped age profiles for the bootMatrix function.
}
\description{
Bootstrapping the age profiles of a given PhyloExpressionSet or DivergenceExpressionSet allows
to compute a statistical confidence that an observed phylotranscriptomic pattern is significant or not.
This function fastly computes for a given PhyloExpressionSet or DivergenceExpressionSet n age profiles with
permuted (random) age assignment of genes.

}
\usage{
cpp_bootMatrix(ExpressionMatrix, AgeVector, permutations)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ExpressionMatrix}{
a numeric matrix specifying a standard expression matrix.
}
  \item{AgeVector}{
a numeric vector storing the age information for each corresponding gene within the expression matrix.
}
  \item{permutations}{
a numeric value specifying the numer of permutations that shall be performed.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a numeric matrix storing n permuted age profiles.
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
\code{\link{bootMatrix}}
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }

