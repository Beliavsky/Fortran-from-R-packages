\encoding{latin1}

\name{HR.Mest}
\alias{HR.Mest}

\title{Simultaneous
Affine Equivariant Estimation of Multivariate Median and Tyler's Shape
Matrix}
\description{
  iterative algorithm that finds the affine equivariant multivariate median by estimating \code{\link{tyler.shape}} simultaneously.
}
\usage{
HR.Mest(X, maxiter = 100, eps.scale = 1e-06, eps.center = 1e-06,
        na.action = na.fail)
}

\arguments{
  \item{X}{a numeric data frame or matrix.  }
  \item{maxiter}{maximum number of iterations. }
  \item{eps.scale}{convergence tolerance for the Tyler's shape matrix subroutine.}
  \item{eps.center}{convergence tolerance for the location estimate.}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
}
\details{
  The algorithm follows the idea of Hettmansperger and Randles (2002). There are, however, some differences. This algorithm
  has the vector of marginal medians as starting point for the location and the starting shape matrix is Tyler's shape
  matrix based on the vector of marginal medians and has then a location step and a shape step which are:
  \describe{
  \item{location step k+1:}{transforming the data as \eqn{y=x V_{k}^{-\frac{1}{2}}}{y = x V_k^-0.5} and
                       computing the spatial median \eqn{\mu_y}{mu_y} of y using the function \code{spatial.median}. Then
                       retransforming \eqn{\mu_y}{mu_y} to the original scale \eqn{\mu_{x,k+1}=\mu_y V_{k}^{\frac{1}{2}} }{mu_x,k+1 = mu_y V_k^0.5 }.
                      }
  \item{shape step k+1:}{computing Tyler's shape matrix  \eqn{V_{k+1}}{V_k+1} with respect to \eqn{\mu_{x,k+1}}{mu_x,k+1} by using the function \code{tyler.shape}.}
  }
  The algorithm stops when the difference between two subsequent location estimates is smaller than \code{eps.center}.

  There is no proof that the algorithm converges.
}
\value{
 A list containing:
  \item{center}{vector with the estimated loaction.}
  \item{scatter}{matrix of the estimated scatter.}
}
\references{\cite{Hettmansperger, T.P. and Randles, R.H. (2002), A practical affine equivariant multivariate median, \emph{Biometrika}, \bold{89},  851--860. }}

\author{Klaus Nordhausen and Seija Sirkia}

\examples{
set.seed(654321)
cov.matrix <- matrix(c(3,2,1,2,4,-0.5,1,-0.5,2), ncol=3)
X <- rmvnorm(100, c(0,0,0), cov.matrix)
res <- HR.Mest(X)
colMeans(X)
res$center
cov.matrix/det(cov.matrix)^(1/3)
res$scatter
rm(.Random.seed)
         }


\keyword{multivariate}
\keyword{nonparametric}
