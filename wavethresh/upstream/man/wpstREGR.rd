\name{wpstREGR}
\alias{wpstREGR}
\title{Construct data frame using new time series using information
	from a previously constructed wpstRO object
}
\description{The \code{\link{makewpstRO}} function takes two time series,
performs a nondecimated wavelet packet transform with the "dependent"
variable one, stores the "best" packets (those that individually correlate
with the response series) and returns the data frame that contains the
response and the best packets. The idea is that the user then performs
some kind of modelling between response and packets. This function takes
a new "dependent" series and returns the best packets in a new data frame
in the same format as the old one. The idea is that the model and the new
data frame can be used together to predict new values for the response
}
\usage{
wpstREGR(newTS, wpstRO)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{newTS}{The new "dependent" time series
}
  \item{wpstRO}{The previously constructed wpstRO object made by
	\code{\link{makewpstRO}}
}
}
\details{
Description says it all
}
\value{
New values of the response time series
}
\references{
See reference to Nason and Sapatinas paper in the help for
\code{\link{makewpstRO}}.
}
\author{
G P Nason
}

\seealso{\code{\link{makewpstRO}}, \code{\link{wpst}}
}
\examples{
#
# See extended example in makewpstRO help, includes example of using this fn
#
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{ts}
