\name{margcond}
\alias{marg_mvnorm}
\alias{cond_mvnorm}
\title{
    Marginal and Conditional Multivariate Normal Distributions
}
\description{
    Computes means and Cholesky factors of covariance or precision matrices
    of multiple multivariate normal distributions.
}
\usage{
marg_mvnorm(chol, invchol, which = 1L)
cond_mvnorm(chol, invchol, which_given = 1L, given, center = FALSE)
}
\arguments{
  \item{chol}{Cholesky factors of covariance matrices as
    \code{\link{ltMatrices}} object, length is recylced to length of \code{lower} and
\code{upper}.
}
  \item{invchol}{Cholesky factors of precision matrices as
    \code{\link{ltMatrices}} object, length is recylced to length of \code{lower} and
\code{upper}. Either \code{chol} or \code{invchol} must be given.
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
  \item{center}{logical, if \code{TRUE}, the negative rescaled conditional mean is returned
   (such that it can be specified as \code{center} argument to
   \code{\link{slpmvnorm}}). By default, the conditional mean is returned.
}
}
\details{
  Derives parameters of the requested marginal or conditional distributions,
  defined by \code{chol} (Cholesky factor of covariance) or \code{invchol}
  (Cholesky factor of precision matrix) and, for conditional distributions,
  the mean.

  More details can be found in the \code{lmvnorm_src} package vignette.

}
\value{
A named list.
}
\seealso{\code{vignette("lmvnorm_src", package = "mvtnorm")}}
\keyword{distribution}
