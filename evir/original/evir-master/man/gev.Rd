\name{gev}
\alias{gev}
\title{Fit Generalized Extreme Value Distribution}
\description{
  Fits generalized extreme value distribution (GEV) to
  block maxima data. 
}
\usage{
gev(data, block = NA, \dots)
}
\arguments{
  \item{data}{data vector. Interpretation depends on value of
    block: if no block size is specified then data are
    interpreted as block maxima; if block size is set, then data
    are interpreted as raw data and block maxima are calculated.}
  \item{block}{the block size. A numeric value is interpreted as the
    number of data values in each successive block. All the data is
    used, so the last block may not contain \code{block} observations.
    If the \code{data} has a \code{times} attribute containing (in
    an object of class \code{"POSIXct"}, or an object that can be
    converted to that class; see \code{\link{as.POSIXct}}) the
    times/dates of each observation, then \code{block} may instead
    take the character values \code{"month"}, \code{"quarter"},
    \code{"semester"} or \code{"year"}.}
  \item{\dots}{arguments passed to \code{optim}}
}
\value{
  An object of class \code{gev} describing the fit and including
  parameter estimates and standard errors.
  Fitting is carried out using maximum likelihood.
}
\seealso{\code{\link{plot.gev}}, \code{\link{gumbel}},
  \code{\link{optim}}, \code{\link{as.POSIXct}}}

\examples{
# Fit GEV to monthly maxima
data(bmw)
out <- gev(bmw, "month") 
# Fit GEV to maxima of blocks of 100 observations
out <- gev(bmw, 100) 
# Fit GEV to the data in nidd.annual, the annual maximum water 
# levels of the River Nidd, using the "BFGS" optimization method
data(nidd.annual) 
out <- gev(nidd.annual, method = "BFGS", control = list(maxit = 500))   
}
\keyword{models}
