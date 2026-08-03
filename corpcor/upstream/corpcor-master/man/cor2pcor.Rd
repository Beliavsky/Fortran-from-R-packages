\name{cor2pcor}
\alias{cor2pcor}
\alias{pcor2cor}

\title{Compute Partial Correlation from Correlation Matrix -- and Vice Versa}
\usage{
cor2pcor(m, tol)
pcor2cor(m, tol)
}
\arguments{
  \item{m}{covariance matrix or (partial) correlation matrix}   
  \item{tol}{tolerance - singular values larger than
             tol are considered non-zero (default value:
	     \code{tol = max(dim(m))*max(D)*.Machine$double.eps}).
	     This parameter is needed for  the singular
	     value decomposition on which \code{\link{pseudoinverse}} is based.
	     }
}
\description{
  \code{cor2pcor} computes the pairwise 
  \emph{partial} correlation coefficients from either a correlation 
  or a covariance matrix. 

  \code{pcor2cor} takes either a partial correlation matrix or 
   a partial covariance matrix as input,
   and computes from it the corresponding correlation matrix.
}

\details{
  The partial
  correlations are the negative standardized concentrations (which in 
  turn are the off-diagonal elements of the inverse correlation or 
  covariance matrix).  In graphical Gaussian models the partial 
  correlations represent the 
  direct interactions between two variables, conditioned on all
  remaining variables.
  
  In the above functions the \code{\link{pseudoinverse}} is employed
  for inversion - hence even singular covariances (with some
  zero eigenvalues) may be used.  However, a better option may be to
  estimate a positive definite covariance matrix using 
  \code{\link{cov.shrink}}.
  
  Note that for efficient computation of partial correlation coefficients from 
  data  \code{x} it is advised to use  \code{pcor.shrink(x)}  and \emph{not} 
  \code{cor2pcor(cor.shrink(x))}.
}

\value{
  A matrix with the pairwise partial correlation coefficients
  (\code{cor2pcor}) or with pairwise
  correlations (\code{pcor2cor}).
}


\author{
  Korbinian Strimmer (\url{https://strimmerlab.github.io}).
}
\references{
  Whittaker J. 1990.  Graphical Models in Applied Multivariate Statistics.
   John Wiley, Chichester.
}

\seealso{\code{\link{decompose.invcov}}, \code{\link{pcor.shrink}}, \code{\link{pseudoinverse}}.}

\examples{
# load corpcor library
library("corpcor")

# covariance matrix
m.cov = rbind(
 c(3,1,1,0),
 c(1,3,0,1),
 c(1,0,2,0),
 c(0,1,0,2)
)
m.cov


# corresponding correlation matrix
m.cor.1 = cov2cor(m.cov)
m.cor.1

# compute partial correlations (from covariance matrix)
m.pcor.1 = cor2pcor(m.cov)
m.pcor.1

# compute partial correlations (from correlation matrix)
m.pcor.2 = cor2pcor(m.cor.1)
m.pcor.2



zapsmall( m.pcor.1 ) == zapsmall( m.pcor.2 )

# backtransformation
m.cor.2 = pcor2cor(m.pcor.1)
m.cor.2
zapsmall( m.cor.1 ) == zapsmall( m.cor.2 )

}
\keyword{multivariate}
