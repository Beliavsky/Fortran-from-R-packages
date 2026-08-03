\name{shrink.intensity}
\alias{estimate.lambda}
\alias{estimate.lambda.var}


\title{Estimation of Shrinkage Intensities}
\description{
  The functions \code{estimate.lambda} and \code{estimate.lambda.var}  
  shrinkage intensities used for correlations and variances used
  in \code{\link{cor.shrink}} and \code{\link{var.shrink}}, respectively.
}

\usage{
estimate.lambda(x, w, verbose=TRUE)
estimate.lambda.var(x, w, verbose=TRUE)
}

\arguments{
  \item{x}{a data matrix}
  \item{w}{optional: weights for each data point - if not specified uniform weights are assumed
            (\code{w = rep(1/n, n)} with \code{n = nrow(x)}).}
  \item{verbose}{report shrinkage intensities (default: TRUE)}
}
\details{
   
   \code{var.shrink} computes the empirical variance of each considered random variable, 
   and shrinks them towards their median.  The corresponding 
   shrinkage intensity \code{lambda.var} is estimated using
    
      \deqn{\lambda_{var}^{*} = (  \sum_{k=1}^p Var(s_{kk}) )/ \sum_{k=1}^p (s_{kk} - median(s))^2  }

  where \eqn{median(s)} denotes the median of the empirical variances (see Opgen-Rhein and Strimmer 2007).   
  
  Similarly,   \code{cor.shrink} computes a shrinkage estimate of the correlation matrix by
  shrinking the empirical correlations towards  the identity matrix.  
  In this case the shrinkage intensity \code{lambda} equals
   
      \deqn{\lambda^{*} = \sum_{k \neq l} Var(r_{kl}) / \sum_{k \neq l} r_{kl}^2  }
  (Sch\"afer and Strimmer 2005).

  Ahdesm\"aki suggested (2012) a computationally highly efficient algorithm to compute 
  the shrinkage intensity estimate for the correlation matrix (see the R code for the implementation).
 
}
\value{
  \code{estimate.lambda} and \code{estimate.lambda.var} returns a number between 0 and 1. 
}
\author{
  Juliane Sch\"afer,
  Rainer Opgen-Rhein, Miika Ahdesm\"aki
  and Korbinian Strimmer (\url{https://strimmerlab.github.io}).
  
}
\references{

  Opgen-Rhein, R., and K. Strimmer. 2007. Accurate ranking of 
  differentially expressed genes by a distribution-free shrinkage 
  approach.    Statist. Appl. Genet. Mol. Biol. \bold{6}:9.
  <DOI:10.2202/1544-6115.1252>

  Sch\"afer, J., and K. Strimmer. 2005.  A shrinkage approach to large-scale
  covariance estimation and implications for functional genomics. 
  Statist. Appl. Genet. Mol. Biol. \bold{4}:32.
  <DOI:10.2202/1544-6115.1175>

}
\seealso{\code{\link{cor.shrink}}, \code{\link{var.shrink}}.}

\examples{
# load corpcor library
library("corpcor")

# small n, large p
p = 100
n = 20

# generate random pxp covariance matrix
sigma = matrix(rnorm(p*p),ncol=p)
sigma = crossprod(sigma)+ diag(rep(0.1, p))

# simulate multinormal data of sample size n  
sigsvd = svd(sigma)
Y = t(sigsvd$v \%*\% (t(sigsvd$u) * sqrt(sigsvd$d)))
X = matrix(rnorm(n * ncol(sigma)), nrow = n) \%*\% Y


# correlation shrinkage intensity
estimate.lambda(X) 
c = cor.shrink(X)
attr(c, "lambda")

# variance shrinkage intensity
estimate.lambda.var(X) 
v = var.shrink(X)
attr(v, "lambda.var")

}
\keyword{multivariate}
