\name{std.error}
\alias{std.error}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to compute the standard error of a numeric vector.
}
\description{
This function computes the standard error of a given numeric vector.}
\usage{
std.error(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a numeric vector for which a standard error should be computed.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a numeric value representing the standard error of a given numeric vector.
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
\code{\link{cpp_std_error}}
}
\examples{

# a simple numeric vector
x <- 1:20

# compute the standard error of x
std.error(x)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }

