\name{print.w2m}
\alias{print.w2m}
\title{
Print a w2m class object
}
\description{
These objects are the matrix representation of a nondecimated wavelet
packet object
}
\usage{
\method{print}{w2m}(x, maxbasis = 10, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	The w2m object to print
}
  \item{maxbasis}{
	The maximum number of basis functions to report on
}
  \item{\dots}{
	Other arguments
}
}
\details{
Prints out information about a w2m object. This function gets called
during \code{\link{makewpstRO}}, and so you can see its output in the
example code in that help function
}
\value{
None
}
\author{
G P Nason
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{makewpstRO}},\code{\link{wpst2m}}}
\examples{
#
# See example in makewpstRO
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{print}
