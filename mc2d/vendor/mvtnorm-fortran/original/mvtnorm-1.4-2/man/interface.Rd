\name{interface}
\alias{mvnorm}
\alias{aperm.mvnorm}
\alias{simulate.mvnorm}
\alias{logLik.mvnorm}
\alias{lLgrad}
\alias{lLgrad.mvnorm}
\alias{margDist}
\alias{margDist.mvnorm}
\alias{condDist}
\alias{condDist.mvnorm}
\title{
    User Interface to Multiple Multivariate Normal Distributions
}
\description{
    A simple user interface for computing on multiple multivariate
    normal distributions.
}
\usage{
mvnorm(mean, invcholmean, chol, invchol)
\S3method{aperm}{mvnorm}(a, perm, ...)
margDist(object, which, ...)
\S3method{margDist}{mvnorm}(object, which, ...)
condDist(object, which_given, given, ...)
\S3method{condDist}{mvnorm}(object, which_given, given, ...)
\S3method{simulate}{mvnorm}(object, nsim = dim(object$scale)[1L], seed = NULL, 
                            standardize = FALSE, as.data.frame = FALSE, ...)
\S3method{logLik}{mvnorm}(object, obs, lower, upper, standardize = FALSE, ...)
\S3method{lLgrad}{mvnorm}(object, obs, lower, upper, standardize = FALSE, ...)
}
\arguments{
  \item{chol}{either an \code{ltMatrices} object specifying (multiple)
              Cholesky factors of the covariance matrix or
              one single numeric lower triangular square matrix.
}
  \item{invchol}{either an \code{ltMatrices} object specifying (multiple)
              inverse Cholesky factors of the covariance matrix or
              one single numeric lower triangular square matrix.
}
  \item{a,object}{objects of class \code{mvnorm}.
}
  \item{perm}{a permutation of the covariance matrix corresponding to \code{a}.
}
  \item{which}{names or indices of elements those marginal distribution
               is of interest.
}
  \item{which_given}{names or indices of elements to condition on.
}
\item{given}{matrix of realisations to condition on (number of rows is
             equal to \code{length(which)}, the number of 
             columns corresponds to the number of matrices in \code{chol}
             or \code{invchol}.
}
 \item{lower}{matrix of lower limits (one column for each observation, \eqn{J} rows).
}
  \item{upper}{matrix of upper limits (one column for each observation, \eqn{J} rows).
}
  \item{obs}{matrix of exact observations (one column for each observation, \eqn{J} rows).
}
  \item{mean}{matrix of means (one column for each observation, length is
             recycled to length of \code{obs}, \code{lower} and \code{upper}).
}
  \item{invcholmean}{matrix of scaled means, that is, \code{invchol \%*\% mean}, one column for each observation, length is
             recycled to length of \code{obs}, \code{lower} and \code{upper}.
}
\item{seed}{an object specifying if and how the random number generator
            should be initialized, see \code{\link[stats]{simulate}}. 
}
\item{standardize}{logical, should the Cholesky factor (or its inverse) undergo
                   standardization (ensuring the covariance matrix is a correlation
                   matrix) before computing the likelihood.
}
\item{nsim}{number of samples to draw.
}
\item{as.data.frame}{logical, convert the $J x N$ matrix result to a
                     classical $N x J$ data frame.
}
\item{\dots}{Additional arguments to \code{\link{ldpmvnorm}} and
             \code{\link{sldpmvnorm}}
}
}
\details{
    The constructor \code{mvnorm} can be used to specify (multiple)
    multivariate normal distributions. \code{margDist} derives marginal and
    \code{condDist} conditional distributions from such objects. A
    \code{simulate} method exists for drawn samples from multivariate
    normals.

    The continuous (data in \code{obs}), discrete (intervals in \code{lower}
    and \code{upper}), and mixed continuous-discrete log-likelihood is
    implemented in \code{logLik}. The corresponding gradients with respect
    to all model parameters and with respect to the data arguments
    is available from \code{lLgrad}.

    Rationals and examples are given in Chapter 7 of the package vignette
    linked to below.
}
\value{
    \code{mvnorm}, \code{margDist}, and \code{condDist} return objects
    of class \code{mvnorm}. \code{logLik} returns the log-likelihood 
    and \code{lLgrad} a list with gradients.
}
\seealso{\code{vignette("lmvnorm_src", package = "mvtnorm")}}
\keyword{distribution}
