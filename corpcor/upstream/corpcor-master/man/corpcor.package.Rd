\encoding{latin1}
\name{corpcor-package}
\alias{corpcor-package}
\docType{package}
\title{The corpcor Package}

\description{
This package implements a James-Stein-type shrinkage estimator for 
  the covariance matrix, with separate shrinkage for variances and correlations.  
  The details of the method are explained in Sch\"afer and Strimmer (2005) 
  <DOI:10.2202/1544-6115.1175> and Opgen-Rhein and Strimmer (2007) 
  <DOI:10.2202/1544-6115.1252>.  The approach is both computationally as well
  as statistically very efficient, it is applicable to ``small n, large p'' data, 
  and always returns a positive definite and well-conditioned covariance matrix.  
  In addition to inferring the covariance matrix the package also provides 
  shrinkage estimators for partial correlations, partial variances, and 
  regression coefficients.  The inverse of the covariance and correlation matrix 
  can be efficiently computed, and as well as any arbitrary power of the 
  shrinkage correlation matrix.  Furthermore, functions are available for fast 
  singular value decomposition, for computing the pseudoinverse, and for 
  checking the rank and positive definiteness of a matrix.

  The name of the package refers to \bold{cor}relations and 
  \bold{p}artial \bold{cor}relations.
}

\author{Juliane Sch\"afer, Rainer Opgen-Rhein, Verena Zuber, Miika Ahdesm\"aki, 
 A. Pedro Duarte Silva, and Korbinian Strimmer (\url{https://strimmerlab.github.io/})}
\references{
  See website: \url{https://strimmerlab.github.io/software/corpcor/}
}
\keyword{multivariate}
\seealso{
\code{\link{cov.shrink}, \link{invcov.shrink}, \link{powcor.shrink}, \link{pcor.shrink},
 \link{fast.svd}.}
}
