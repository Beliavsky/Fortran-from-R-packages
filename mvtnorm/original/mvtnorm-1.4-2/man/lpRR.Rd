\name{lpRR}
\alias{lpRR}
\alias{slpRR}
\title{
    Multivariate Normal Log-likelihood and Score Functions for Reduced Rank
    Covariances
}
\description{
    Computes the log-likelihood (contributions) of  
    interval-censored observations from multivariate 
    normal distributions with reduced rank structure 
    and evaluates corresponding score functions.
}
\usage{
lpRR(lower, upper, mean = 0, B, D = rep(1, nrow(B)), 
     Z, weights = 1 / ncol(Z), log.p = TRUE)
slpRR(lower, upper, mean = 0, B, D = rep(1, nrow(B)), 
      Z, weights = 1 / ncol(Z), log.p = TRUE)
}
\arguments{
  \item{lower}{vector of lower limits (one element for each dimension, \eqn{J}
elements).
}
  \item{upper}{vector of upper limits (one element for each dimension, \eqn{J}
elements).
}
  \item{mean}{vector of means (one element for each dimension, length is
recycled to length of \code{lower} and \code{upper}).
}
  \item{B}{matrix of dimension \eqn{J \times K}.
}
  \item{D}{vector of \eqn{J} diagonal elements.
}
  \item{Z}{matrix of standard normal random variables, with \eqn{K} nrows.
}
  \item{weights}{optional weights.
}
  \item{log.p}{logical. By default, log-probabilities are returned.
}
}
\details{
  Evaluates the multivariate normal log-likelihood defined by \code{mean},
  \code{B} and \code{D} when the covariance is \eqn{\Sigma = B B^\top + D} 
  over boxes defined by \code{lower} and \code{upper}. Details are given
  in \bibcitet{mvtnorm::Genz_Bretz_2009}, Chapter 2.3.1.
  
  \code{slpRR} computes 
  the corresponding score functions with respect to \code{lower},
  \code{upper}, \code{mean}, \code{B} and \code{D}.

  More details can be found in the \code{lmvnorm_src} package vignette.

}
\value{
The log-likelihood (\code{log.p = TRUE}) or corresponding probability.
\code{slpRR} return the scores.
}
\seealso{\code{vignette("lmvnorm_src", package = "mvtnorm")}}
\references{\bibshow{*}}
\examples{
  J <- 6
  K <- 3
  B <- matrix(rnorm(J * K), nrow = J)
  D <- runif(J)
  S <- tcrossprod(B) + diag(D)
  a <- -(2 + runif(J))
  b <- 2 + runif(J)
  M <- 1e4
  Z <- matrix(rnorm(K * M), nrow = K)
  ## log-likelihood
  lpRR(lower = a, upper = b, B = B, D = D, Z = Z)
  ## score wrt all arguments
  slpRR(lower = a, upper = b, B = B, D = D, Z = Z)
}
\keyword{distribution}
