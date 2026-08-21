\name{MLE for Poisson lognormal distribution}
\alias{poilogMLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Maximum Likelihood Estimation for Poisson Lognormal Distribution }
\description{
  \code{poilogMLE} fits the Poisson lognormal distribution to data and estimates 
  parameters mean \code{mu} and standard deviation \code{sig} in the lognormal distribution
}
\usage{
poilogMLE(n, startVals = c(mu=1, sig=2), 
          nboot = 0, zTrunc = TRUE,
          method = "BFGS", control = list(maxit=1000))
}
\arguments{
  \item{n}{ A vector of counts }
  \item{startVals}{ Starting values of parameters, see details }
  \item{nboot}{ Number of parametric bootstraps, defaults to zero }
  \item{zTrunc}{ Zero-truncated distribution, defaults to TRUE }
  \item{method}{ Method to use during optimization, see details }
  \item{control}{ A list of control parameters for the optimization routine, see details }
  
}
\details{
 
  The function estimates parameters mean \code{mu} and standard deviation \code{sig}. 
  In cases of incomplete sampling the estimate of \code{mu} will be confounded with the sampling
  intensity (see \code{\link{rpoilog}}). Assuming sampling intensity \eqn{\nu}, 
  the estimates of the mean is \eqn{\code{mu}+\ln(\nu)}. Parameter \code{sig} can be estimated without 
  any knowledge of sampling intensity.\cr
  The parameters must be given starting values for the optimization procedure (default starting values are 
  used if starting values are not specified in the function call).
  
  The function uses the optimization procedures in \code{\link{optim}} to obtain the maximum likelihood estimate.
  The \code{method} and \code{control} arguments are passed to \code{\link{optim}}, see the help page for this 
  function for additional methods and control parameters.
    
  A zero-truncated distribution (see \code{\link{dpoilog}}) is assumed by default (\code{zTrunc = TRUE}).
  In cases where the number of zeros is known the \code{zTrunc} argument should be set to \code{FALSE}.
  
  The approximate fraction of species revealed by the sample is \eqn{1-q(0;\code{mu},\code{sig})}.
  
  Parametric bootstrapping is done by simulating new sets of observations using the  estimated parameters 
  (see \code{\link{rbipoilog}}).
  
  
}
\value{
  \item{par }{Maximum likelihood estimates of the parameters}
  \item{p }{Approximate fraction of species revealed by the sample}
  \item{logLval }{Log likelihood of the data given the estimated parameters}
  \item{gof }{Goodness of fit measure obtained by checking the rank of logLval against logLval's obtained 
                during the bootstrap procedure, (gof<0.05) or (gof>0.95) indicates lack of fit}
  \item{boot }{A data frame containing the bootstrap replicates of parameters and logLval}
}
\references{
Bulmer, M. G. 1974. On fitting the Poisson lognormal distribution to species abundance data. Biometrics \bold{30}, 651-660.\cr
Engen, S., R. Lande, T. Walla and P. J. DeVries. 2002. Analyzing spatial structure of communities using the two-dimensional 
            Poisson lognormal species abundance model. American Naturalist \bold{160}, 60-73.}
\author{ Vidar Grotan \email{vidar.grotan@ntnu.no}, Steinar Engen}

\seealso{ \code{\link{optim}}, \code{\link{dpoilog}}, \code{\link{rpoilog}}  }
\examples{
### simulate observations
\dontshow{set.seed(200)}
n <- rpoilog(S=80,mu=1,sig=2)

### obtain estimates of parameters
est <- poilogMLE(n)

### similar, but now with bootstrapping ###
\dontrun{est <- poilogMLE(n,nboot=10)}

### change start values and request tracing information 
### from optimization procedure
est <- poilogMLE(n,startVals=c(2,3),
                 control=list(maxit=1000,trace=1, REPORT=1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
