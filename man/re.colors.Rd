\name{re.colors}
\alias{re.colors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to get a vector of length n storing a palette of colors.
}
\description{
This function returns distinct color names that can be used to color different curves,
so that they can be distiguished quite clearly.}
\usage{
re.colors(n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
the number of colors to be in the palette.
}
}
\details{
}
\value{
a character vector containing different color names that can be used for relative expression plots. 
}
\references{
\code{\link{palette}}
}
\author{
Hajk-Georg Drost
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{bar.colors}}
}
\examples{

# get 12 different colors for 12 different phylostrata
re_colors <- re.colors(12)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }

