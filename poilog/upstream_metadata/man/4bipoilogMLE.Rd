\name{MLE for bivariate Poisson lognormal distribution}
\alias{bipoilogMLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Maximum Likelihood Estimation for Bivariate Poisson Lognormal Distribution }
\description{
  \code{bipoilogMLE} fits the bivariate Poisson lognormal distribution to data
}
\usage{
bipoilogMLE(n1, n2 = NULL, 
            startVals = c(mu1=1, mu2=1, sig1=2, sig2=2, rho=0.5), 
            nboot = 0, zTrunc = TRUE, file = NULL, 
            method = "BFGS", control = list(maxit=1000))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n1}{ a vector or a matrix with two columns of pairwise counts of observed individuals for each species }
  \item{n2}{ if n1 is not given as a matrix, a vector of counts with same ordering of species as in argument n1}
  \item{startVals}{ starting values of parameters }
  \item{nboot}{ number of parametric bootstraps, defaults to zero }
  \item{zTrunc}{ logical; if TRUE (default) the zero-truncated distribution is fitted }
  \item{file}{ text file to hold copies of bootstrap estimates }
  \item{method}{ method to use during optimization, see details }
  \item{control}{ a list of control parameters for the optimization routine, see details }
}
\details{
  
  The function estimates the parameters \code{mu1}, \code{sig1}, \code{mu2}, \code{sig2} and \code{rho}. 
  In cases of incomplete sampling the estimates of \code{mu1} and \code{mu2} will be confounded with the sampling
  intensities (see \code{\link{rbipoilog}}). Assuming sampling intensities \eqn{\nu_1}{\nu1} and \eqn{\nu_2}{\nu2}, 
  the estimates of the means are \eqn{\code{mu1}+\ln \nu_1}{\code{mu1}+ \ln nu1} and \eqn{\code{mu2}+\ln\nu_2}{\code{mu2}+ \ln nu2}. Parameters
  \code{sig1}, \code{sig2} and \code{rho} can be estimated without any knowledge of sampling intensities.
  The parameters must be given starting values for the optimization procedure (default starting values are 
  used if starting values are not specified in the function call).
  
  A zero-truncated distribution (see \code{\link{dbipoilog}}) is assumed by default (\code{zTrunc = TRUE}).
  In cases where the number of zeros is known the \code{zTrunc} argument should be set to \code{FALSE}.
  
  The function uses the optimization procedures in \code{\link{optim}} to obtain the maximum likelihood estimate.
  The \code{method} and \code{control} arguments are passed to \code{\link{optim}}, see the help page for this 
  function for additional methods and control parameters.
  
  The approximate fraction of species revealed by each sample is estimated as 1 minus the zero term of the 
  univariate Poisson lognormal distribution: \eqn{1-q(0;\code{mu1},\code{sig1})} and \eqn{1-q(0;\code{mu2},\code{sig2})}.

  Parametric bootstrapping could be time consuming for large data sets. If argument \code{file} is specified, e.g. 
  \code{file =}\sQuote{\code{C:\\\\myboots.txt}}, the matrix with bootstrap estimates are copied into a tab-separated text-file 
  providing extra backup. Bootstrapping is done by simulating new sets of observations conditioned on the observed number of 
  species (see \code{\link{rbipoilog}}).
  
}
\value{
  \item{par }{Maximum likelihood estimates of the parameters}
  \item{p }{Approximate fraction of species revealed by the samples for sample 1 and 2 respectively}
  \item{logLval }{Log likelihood of the data given the estimated parameters}
  \item{gof }{Goodness of fit measure obtained by checking the rank of logLval against logLval's obtained 
              during the bootstrap procedure, (\code{gof} < 0.05) or (\code{gof} > 0.95) indicates lack of fit}
  \item{boot }{A data frame containing the bootstrap replicates of parameters and logLval}
  
}
\references{ 
Engen, S., R. Lande, T. Walla and P. J. DeVries. 2002. Analyzing spatial structure of communities using the two-dimensional 
            Poisson lognormal species abundance model. American Naturalist \bold{160}, 60-73.}
\author{ Vidar Grotan \email{vidar.grotan@ntnu.no}, Steinar Engen   }

\seealso{ \code{\link{optim}}, \code{\link{dbipoilog}}, \code{\link{rbipoilog}}  }
\examples{
\dontshow{set.seed(200)}
## simulate observations
xy  <- rbipoilog(S=30,mu1=1,mu2=1,sig1=2,sig2=2,rho=0.5)

## obtain estimates of parameters
est <- bipoilogMLE(xy)

## change start values and request tracing information 
## from optimization procedure
\dontrun{est <- bipoilogMLE(xy,startVals=c(2,2,4,4,0.3),
       control=list(maxit=1000,trace=1, REPORT=1))}

## effect of sampling intensity 
xy  <- rbipoilog(S=100,mu1=1,mu2=1,sig1=2,sig2=2,rho=0.5,nu1=0.5,nu2=0.5)
\dontrun{est <- bipoilogMLE(xy)}
## the expected estimates of mu1 and mu2 are now 1-log(0.5) = 0.3 (approximately) 


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
