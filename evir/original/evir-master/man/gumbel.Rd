\name{gumbel}
\alias{gumbel}
\title{Fit Gumbel Distribution}
\description{
  Fits gumbel distribution (GEV with \code{xi = 0}) to
  block maxima data. 
}
\usage{
gumbel(data, block = NA, \dots)
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
  An object of class \code{"gev"} describing the fit and including
  parameter estimates and standard errors.
  Fitting is carried out using maximum likelihood.
}
\details{ 
  This function is primarily intended for comparison with GEV for 
  assessing the need for a heavy-tailed Frechet (or short-tailed 
  Weibull) to model block maxima. 
}
\seealso{\code{\link{plot.gev}}, \code{\link{gev}},
  \code{\link{optim}}, \code{\link{as.POSIXct}}}

\examples{
# Fit Gumbel to maxima of blocks of 100 observations
data(bmw)
out <- gumbel(bmw, 100) 
# Fit Gumbel to the data in nidd.annual, the annual maximum water 
# levels of the River Nidd, using the "BFGS" optimization method
data(nidd.annual)
out <- gumbel(nidd.annual, method = "BFGS", control = list(maxit = 500))
}
\keyword{models}
