\name{pot}
\alias{pot}
\title{Peaks Over Thresholds Model}
\description{
  Fits a Poisson point process to the data, an approach 
  sometimes known as peaks over thresholds (POT), and
  returns an object of class \code{"potd"}. 
}
\usage{
pot(data, threshold = NA, nextremes = NA, run = NA, picture = TRUE,
    \dots)
}
\arguments{
  \item{data}{numeric vector of data, which may have a \code{times}
    attribute containing (in an object of class \code{"POSIXct"}, or
    an object that can be converted to that class; see
    \code{\link{as.POSIXct}}) the times/dates of each observation.
    If no \code{times} attribute exists, the data are assumed to
    be equally spaced.}
  \item{threshold}{a threshold value (either this or \code{nextremes}
    must be given but not both)}
  \item{nextremes}{the number of upper extremes to be used (either
    this or \code{threshold} must be given but not both)}
  \item{run}{if the data are to be declustered the run length
    parameter for the runs method (see \code{\link{decluster}})
    should be entered here}
  \item{picture}{whether or not a picture should be drawn if
    declustering is performed}
  \item{\dots}{arguments passed to \code{optim}}
}
\value{
  An object of class \code{"potd"} describing the fit and including
  parameter estimates and standard errors. 
}
\details{
  Uses \code{optim} for point process likelihood maximization. 
}
\seealso{\code{\link{gpd}}, \code{\link{plot.potd}}, 
  \code{\link{plot.gpd}}, \code{\link{decluster}},
  \code{\link{optim}}, \code{\link{as.POSIXct}}}

\examples{
data(danish)
out <- pot(danish, 10)  
# Fits POT model to Danish fire insurance losses 
}
\keyword{models}
