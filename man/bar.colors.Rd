\name{bar.colors}
\alias{bar.colors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function to get a vector of length n storing a palette of colors for multiple bars in barplots.
}
\description{
A nice color palette for barplots with several bars.
}
\usage{
bar.colors(n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
the number of colors to be in the palette.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a character vector containing different color names that can be used for barplots. 
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
\code{\link{re.colors}}
}
\examples{

# get 5 different colors for 5 different bars
barplot_colors <- bar.colors(5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }

