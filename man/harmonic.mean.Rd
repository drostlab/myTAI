\name{harmonic.mean}
\alias{harmonic.mean}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the harmonic mean of a vector.}
\description{
This function computes the harmonic mean of a given numeric vector x.
}
\usage{
harmonic.mean(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a numeric vector for which the harmonic mean shall be computed.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a numeric value representing the harmonic mean numeric input vector.

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
\code{\link{geom.mean}}
}
\examples{
# a simple vector
x <- 1:20

# the harmonic mean of the simple vector x
harmonic.mean(x)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }

