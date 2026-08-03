\name{HP1.shape}
\alias{HP1.shape}

\title{ One Step Rank Scatter Estimator  }
\description{
  one step M-estimator of the scatter matrix based on ranks.
}
\usage{
HP1.shape(X, location = "Estimate", na.action = na.fail, ...)
}

\arguments{
  \item{X}{a numeric data frame or matrix. }
  \item{location}{ if 'Estimate' the location and scatter matrix used for computing the spatial signs are estimated simultaneously using \code{\link{HR.Mest}},
  if 'Origin' or numeric \code{\link{tyler.shape}} is used with respect to origin or the given value, respectively, to obtain the spatial signs.   }
  \item{na.action}{ a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
  \item{\dots}{ arguments that can be passed on to \code{\link{tyler.shape}} or \code{\link{HR.Mest}}.}
}
\details{
  This is a one step M-estimator of shape which is standardized in such a way that the determinant is 1.

  The exact formula is:
  \deqn{V = V_{0}^{\frac{1}{2}} ave\{a(\frac{R_{i}}{n+1})u_{i}'u_{i} \} V_{0}^{\frac{1}{2}}.}{%
  V=V_0^1/2 ave\{a(R_i / (n+1)) u_i' u_i\} V_0^1/2.}

  where \eqn{V_{0}}{V_0} is Tyler's shape matrix, \eqn{u_{i}=||z_{i}||^{-1} z_{i}}{u_i=||z_{i}||^-1 z_i} is the spatial sign of \eqn{z_{i}=(x_{i}-\mu) V_{0}^{-\frac{1}{2}}}{z_i=(x_i-mu) V_0^-(1/2)}
  and \eqn{R_{i}}{R_i}
  gives the rank of \eqn{||z_{i}||}{||z_i||} among \eqn{||z_{1}||,\ldots,||z_{n}||}{||z_1||,...,||z_n||}. The van der Warden score function \eqn{a(.)}{a(.)} is the inverse of the cdf of a chi-squared distribution with p degrees of freedom.

  This scatter matrix is based on the test for shape developed in the paper by Hallin and Paindaveine (2006), its usage with respect to the origin is demonstrated in Nordhausen et al. (2006).
}

\references{\cite{Hallin, M. and Paindaveine, D. (2006), Semiparametrically efficient rank-based inference for shape. I. Optimal rank-based tests for sphericity, \emph{Annals of Statistics}, \bold{34},  2707--2756. }

            \cite{Nordhausen, K., Oja, H. and Paindaveine, D. (2009), Signed-rank tests for location in the symmetric independent component model, \emph{Journal of Multivariate Analysis}, \bold{100},  821--834. }}


\author{Klaus Nordhausen}

\examples{
set.seed(654321)
cov.matrix <- matrix(c(3,2,1,2,4,-0.5,1,-0.5,2), ncol=3)
X <- rmvnorm(100, c(0,0,0), cov.matrix)
HP1.shape(X)
HP1.shape(X, location="Origin")
cov.matrix/det(cov.matrix)^(1/3)
rm(.Random.seed)
         }


\keyword{ multivariate }
\keyword{nonparametric}
