\name{pow}
\alias{pow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute different power values of a given data vector.}
\description{
This function transforms a given numeric vector (data) into a numeric vector of power N.}
\usage{
pow(x, power)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a numeric vector for which each element the power N shall be computed.
}
  \item{power}{
the power to be computed for the numeric vector. E.g. pow(data,2) = data^2.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a numeric value representing the power transformed value.
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

x <- 1:20
# x^2
pow(x, power = 2)

# x^5
pow(x, power = 2)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
