\name{rlevel.gev}
\alias{rlevel.gev}
\title{Calculate Return Levels Based on GEV Fit}
\description{
  Calculates the k-block return level and 95\% confidence interval
  based on a GEV model for block maxima, where k is specified by
  the user. The k-block return level is that level exceeded once
  every k blocks, on average. 
}
\usage{
rlevel.gev(out, k.blocks = 20, add = FALSE, \dots)
}
\arguments{
  \item{out}{an object returned by the function \code{gev}}
  \item{k.blocks}{specifies the particular return level to be
    estimated; default set arbitrarily to 20}
  \item{add}{whether the return level should be added graphically
    to a time series plot; if \code{FALSE} a graph of the profile
    likelihood curve showing the return level and its confidence
    interval is produced}
  \item{\dots}{other graphics parameters}
}
\value{
  Vector containing lower 95\% bound of confidence interval,
  estimated return level and upper 95\% bound. 
}
\details{
  The GEV likelihood is reparameterized in terms of the unknown
  return level and profile likelihood arguments are used to
  construct a confidence interval. 
}
\seealso{\code{\link{gev}}, \code{\link{plot.gev}}}

\examples{
data(bmw)
out <- gev(bmw, "month")
# Fit GEV to monthly maxima of daily returns on BMW share price
\dontrun{rlevel.gev(out, 40)}
# Calculate the 40 month return level
}
\keyword{htest}
