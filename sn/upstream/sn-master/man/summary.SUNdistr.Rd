%  file sn/man/summary.SUNdistr.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013-2021 Adelchi Azzalini
%---------------------
\name{summary.SUNdistr}
\encoding{UTF-8}

\alias{summary.SUNdistr}
\alias{summary,SUNdistr-method}

\title{Summary of a \acronym{SUN} distribution object}

\description{Produce a summary of an object of class \code{"SUNdistr"}}

\usage{\S4method{summary}{SUNdistr}(object, ...)}

\arguments{
  \item{object}{an object of class \code{"SUNdistr"}.}

  \item{...}{optional arguments passed to \code{mom.mtruncnorm} for the
     regulation of its working.qq}
}

\value{An S4-object with the following slots:

  \item{dp}{the parameters of the distrbution, a list}
  \item{name}{the name of the distribution, a character string}
  \item{compNames}{the names of the components, a character vector}
  \item{HcompNames}{the names of the hidden components, a character vector}   
  \item{mean}{the mean value, a vector}
  \item{var.cov}{the variance-covariance matrix}
  \item{gamma1}{the marginal indices of asymmetry, a vector}
  \item{cum3}{the third order cumulants, a three-dimensional array}
  \item{mardia}{the Mardia's measures of multivariate asymmetry and
     skewness, a vector of length two}
}
 
\author{Adelchi Azzalini}

%\note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{makeSUNdistr}} for building a \acronym{SUN} distribution object
 
 methods \code{\link[base]{mean}} and  \code{\link[stats]{vcov}} 
 for computing the mean vector and the variance matrix of 
 \code{\link{SUNdistr-class}} objects
 
 }
 
 
\references{
Arellano-Valle, R. B. and Azzalini, A. (2021).
Some properties of the unified skew-normal distribution.
\emph{Statistical Papers}, 
\doi{https://doi.org/10.1007/s00362-021-01235-2}
and \href{https://arxiv.org/abs/2011.06316}{arXiv:2011.06316}    

Azzalini, A. with the collaboration of Capitanio, A. (2014). 
  \emph{The Skew-Normal and Related Families}. 
  Cambridge University Press, IMS Monographs series.
}


\examples{
Omega <- matrix(c(5, 1, 1, 6), 2, 2)
Delta <- matrix(c(0.30, 0.50, 0.50, 0.85), 2, 2, byrow=TRUE)
Gamma <- matrix(c(1, 0.18, 0.18, 1), 2, 2)
tau <- c(0.4, -0.8)
dp2 <- list(x=c(1, 0), Omega=Omega, Delta=Delta, tau=tau, Gamma=Gamma)
sun2 <-  makeSUNdistr(dp=dp2, name="SUN2", compNames=c("u", "v"))
s <- summary(sun2)
}

\keyword{multivariate}
\keyword{distribution} 
