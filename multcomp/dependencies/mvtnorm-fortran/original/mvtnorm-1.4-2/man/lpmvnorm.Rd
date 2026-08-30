\name{lpmvnorm}
\alias{lpmvnorm}
\alias{slpmvnorm}
\alias{ldmvnorm}
\alias{sldmvnorm}
\alias{ldpmvnorm}
\alias{sldpmvnorm}
\title{
    Multivariate Normal Log-likelihood and Score Functions
}
\description{
    Computes the log-likelihood (contributions) of multiple exact or 
    interval-censored observations (or a mix thereof) from multivariate 
    normal distributions and evaluates corresponding score functions.
}
\usage{
lpmvnorm(lower, upper, mean, invcholmean, center = NULL, chol, invchol, logLik = TRUE, 
         M = NULL, w = NULL, seed = NULL, tol = .Machine$double.eps, fast = FALSE)
slpmvnorm(lower, upper, mean, invcholmean, center = NULL, chol, invchol, logLik = TRUE, 
          M = NULL, w = NULL, seed = NULL, tol = .Machine$double.eps, fast = FALSE)
ldmvnorm(obs, mean, invcholmean, chol, invchol, logLik = TRUE) 
sldmvnorm(obs, mean, invcholmean, chol, invchol, logLik = TRUE) 
ldpmvnorm(obs, lower, upper, mean, invcholmean, chol, invchol, logLik = TRUE, ...) 
sldpmvnorm(obs, lower, upper, mean, invcholmean, chol, invchol, logLik = TRUE, ...) 
}
\arguments{
  \item{lower}{matrix of lower limits (one column for each observation, \eqn{J} rows).
}
  \item{upper}{matrix of upper limits (one column for each observation, \eqn{J} rows).
}
  \item{obs}{matrix of exact observations (one column for each observation, \eqn{J} rows).
}
  \item{mean}{matrix of means (one column for each observation, length is
recycled to length of \code{obs}, \code{lower} and \code{upper}).
}
  \item{invcholmean}{matrix of means left-multiplied with inverse Cholesky factor 
  (\code{invchol \%*\% mean}, one column for each observation, length is
  recycled to length of \code{obs}, \code{lower} and \code{upper}).
}
  \item{center}{matrix of negative rescaled means (one column for each observation, length is
recycled to length of \code{lower} and \code{upper}) as returned by
\code{cond_mvnorm(..., center = TRUE).}. 
}
  \item{chol}{Cholesky factors of covariance matrices as
    \code{\link{ltMatrices}} object, length is recylced to length of
    \code{obs}, \code{lower} and \code{upper}.
}
  \item{invchol}{Cholesky factors of precision matrices as
    \code{\link{ltMatrices}} object, length is recylced to length of \code{lower} and
\code{upper}. Either \code{chol} or \code{invchol} must be given.
}
  \item{logLik}{logical, if \code{TRUE}, the log-likelihood is returned,
otherwise the individual contributions to the sum are returned.
}
  \item{M}{number of iterations, early stopping based on
estimated errors is NOT implemented.
}
  \item{w}{an optional matrix of weights with \eqn{J - 1} rows. This allows to replace the default
Monte-Carlo procedure \bibcitep{mvtnorm::numerical-:1992} with a quasi-Monte-Carlo approach
\bibcitep{mvtnorm::Genz_Bretz_2002}. Note that the same weights for evaluating the
multivariate normal probability are used for all observations when
\code{ncol(w) == M} is specified. If \code{ncol(w) == ncol(lower) * M}, each
likelihood contribution is evaluated on the corresponding sub-matrix.
If \code{w} is \code{NULL}, different uniform numbers are
drawn for each observation.
}
  \item{seed}{an object specifying if and how the random number generator
          should be initialized, see \code{\link[stats]{simulate}}. Only 
          applied when \code{w} is \code{NULL}.
}
  \item{tol}{tolerance limit, values smaller than \code{tol} are interpreted
as zero.
}
  \item{fast}{logical, if \code{TRUE}, a faster but less accurate version 
              of \code{pnorm} is used internally.
}
  \item{\dots}{additional arguments to \code{lpmvnorm}.
}
}
\details{
  Evaluates the multivariate normal log-likelihood defined by \code{means} and
  \code{chol} over boxes defined by \code{lower} and \code{upper} or for
  exact observations \code{obs}. 

  Monte-Carlo \bibcitep{|mvtnorm::numerical-:1992|the default} and quasi-Monte-Carlo 
  \bibcitep{mvtnorm::Genz_Bretz_2002}
  integration is implemented, the latter with weights obtained, for example,
  from packages \pkg{qrng} or \pkg{randtoolbox}. It is the responsibility of
  the user to ensure a meaningful lattice is used. In case of doubt, use
  plain Monte-Carlo (\code{w = NULL}) or \code{\link{pmvnorm}}.

  \code{slpmvnorm} computes both the individual log-likelihood contributions
  and the corresponding score matrix (of dimension \eqn{J \times (J + 1) / 2 \times N}) if
  \code{chol} contains diagonal elements. Otherwise, the dimension is \eqn{J
  \times (J - 1) / 2 \times N}. The scores for exact or mixed exact-interval 
  observations are computed by \code{sldmvnorm} and \code{sldpmvnorm},
  respectively.

  More details can be found in the \code{lmvnorm_src} package vignette.

}
\value{
The log-likelihood (\code{logLik = TRUE}) or the individual contributions to the log-likelihood.
\code{slpmvnorm}, \code{sldmvnorm}, and \code{sldpmvnorm} return the score
matrices and, optionally (\code{logLik = TRUE}), the individual log-likelihood contributions
as well as scores for \code{obs}, \code{lower}, \code{upper}, and
\code{mean}.
}
\seealso{\code{\link{dmvnorm}}, \code{vignette("lmvnorm_src", package = "mvtnorm")}}
\references{\bibshow{*}}
\examples{

  ### five observations
  N <- 5L
  ### dimension
  J <- 4L

  ### lower and upper bounds, ie interval-censoring
  lwr <- matrix(-runif(N * J), nrow = J)
  upr <- matrix(runif(N * J), nrow = J)

  ### Cholesky factor
  (C <- ltMatrices(runif(J * (J + 1) / 2), diag = TRUE))
  ### corresponding covariance matrix
  (S <- as.array(Tcrossprod(C))[,,1])

  ### plain Monte-Carlo (Genz, 1992)
  w <- NULL
  M <- 25000
  ### quasi-Monte-Carlo (Genz & Bretz, 2002, but with different weights)
  if (require("qrng")) w <- t(ghalton(M * N, J - 1))

  ### log-likelihood
  lpmvnorm(lower = lwr, upper = upr, chol = C, w = w, M = M)

  ### compare with pmvnorm
  exp(lpmvnorm(lower = lwr, upper = upr, chol = C, logLik = FALSE, w = w, M = M))
  sapply(1:N, function(i) pmvnorm(lower = lwr[,i], upper = upr[,i], sigma = S))

  ### log-lik contributions and score matrix
  slpmvnorm(lower = lwr, upper = upr, chol = C, w = w, M = M, logLik = TRUE)

}
\keyword{distribution}
