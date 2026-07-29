\name{gpd}
\alias{gpd}
\title{Fit Generalized Pareto Model}
\description{
  Returns an object of class \code{"gpd"} representing the fit of
  a generalized Pareto model to excesses over a high threshold. 
}
\usage{
gpd(data, threshold = NA, nextremes = NA, method = c("ml", "pwm"),
    information = c("observed", "expected"), \dots)
}
\arguments{
  \item{data}{data vector}
  \item{threshold}{a threshold value (either this or
    \code{nextremes} must be given but not both)}
  \item{nextremes}{the number of upper extremes to be used (either
    this or \code{threshold} must be given but not both)}
  \item{method}{whether parameters should be estimated by the maximum
    likelihood method \code{"ml"} or the probability-weighted moments
    method \code{"pwm"}}
  \item{information}{whether standard errors should be calculated with
    \code{"observed"} or \code{"expected"} information. This only applies
    to the maximum likelihood method; for the probability-weighted moments
    method \code{"expected"} information is used if possible}
  \item{\dots}{arguments passed to \code{optim}}
}
\value{
  An object of class \code{"gpd"} describing the fit and including
  parameter estimates and standard errors. 
}
\details{
  The function uses the general purpose optimization function
  \code{optim} when \code{method = "ml"} is chosen.
}
\references{
  Parameter and  quantile  estimation  for  the  generalized
  Pareto   distribution,   JRM   Hosking   and   JR  Wallis,
  \emph{Technometrics} \bold{29}(3), pages 339-349, 1987.
}
\seealso{\code{\link{plot.gpd}}, \code{\link{shape}}, 
  \code{\link{quant}}, \code{\link{optim}}}     

\examples{
data(danish)
out <- gpd(danish, 10) 
# Fits GPD to excess losses over 10 for the Danish
# fire insurance data 
}
\keyword{models}
