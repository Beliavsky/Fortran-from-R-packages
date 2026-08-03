\encoding{latin1}

\name{spatial.median}
\alias{spatial.median}

\title{ Spatial Median }
\description{
  iterative algorithm to compute the spatial median.
}
\usage{ 
spatial.median(X, init = NULL, maxiter = 500, eps = 1e-06, 
               print.it = FALSE, na.action = na.fail)
}

\arguments{
  \item{X}{a numeric data frame or data matrix.  }
  \item{init}{ Starting value for the alogrihtm, if 'NULL', the vector of marginal medians is used. }
  \item{maxiter}{ maximum number of iterations.  }
  \item{eps}{ convergence tolerance. }
  \item{print.it}{logical. If TRUE prints the number of iterations, otherwise not.}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail. }
}
\details{
  Follows the algorithm of Vardi and Zhang.
}
\value{
  vector of the spatial median.
}
\references{\cite{Vardi, Y. and Zhang, C.-H. (1999), The multivariate L1-median and associated data depth, \emph{PNAS}, \bold{97},  1423--1426. }}

\author{Klaus Nordhausen and Seija Sirkia}

\examples{
set.seed(654321)
cov.matrix <- matrix(c(3,2,1,2,4,-0.5,1,-0.5,2), ncol=3)
X <- rmvnorm(100, c(0,0,0), cov.matrix)
spatial.median(X)
rm(.Random.seed)
         }


\keyword{ multivariate }
\keyword{nonparametric}

