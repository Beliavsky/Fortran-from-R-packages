\name{cov.shrink}
\alias{cov.shrink}
\alias{cor.shrink}
\alias{var.shrink}


\title{Shrinkage Estimates of Covariance and Correlation}
\description{
  The functions \code{var.shrink}, \code{cor.shrink}, and \code{cov.shrink} compute 
  shrinkage estimates of variance, correlation, and covariance, respectively.
}

\usage{
var.shrink(x, lambda.var, w, verbose=TRUE)
cor.shrink(x, lambda, w, verbose=TRUE)
cov.shrink(x, lambda, lambda.var, w, verbose=TRUE)
}

\arguments{
  \item{x}{a data matrix}
  \item{lambda}{the correlation shrinkage intensity (range 0-1).
                If \code{lambda} is not specified (the default) it is estimated
		using an analytic formula from Sch\"afer and Strimmer (2005) 
		- see details below.  
		For \code{lambda=0} the empirical correlations are recovered.}
  \item{lambda.var}{the variance shrinkage intensity (range 0-1). 
                If \code{lambda.var} is not specified (the default) it is estimated
	        using an analytic formula from Opgen-Rhein and Strimmer (2007)
                - see details below.  
	        For \code{lambda.var=0} the empirical variances are recovered.}	 
  \item{w}{optional: weights for each data point - if not specified uniform weights are assumed
            (\code{w = rep(1/n, n)} with \code{n = nrow(x)}).}
  \item{verbose}{output some status messages while computing (default: TRUE)}
}
\details{
   
   \code{var.shrink} computes the empirical variance of each considered random variable, 
   and shrinks them towards their median.  The shrinkage intensity is estimated using 
   \code{\link{estimate.lambda.var}} (Opgen-Rhein and Strimmer 2007).
  
  Similarly   \code{cor.shrink} computes a shrinkage estimate of the correlation matrix by
  shrinking the empirical correlations towards  the identity matrix.  
  In this case the shrinkage intensity is computed using  \code{\link{estimate.lambda}}
  (Sch\"afer and Strimmer 2005).
      
   In comparison with the standard empirical estimates 
  (\code{\link{var}}, \code{\link{cov}}, and \code{\link{cor}}) the shrinkage estimates exhibit
   a number of favorable properties.  For instance, 
  \enumerate{
    \item they are typically much more efficient, i.e. they show (sometimes dramatically) better 
          mean squared error,
    \item the estimated covariance and correlation matrices are always  positive definite
          and well conditioned (so that there are no numerical problems when computing their inverse),
    \item  they are inexpensive to compute, and
    \item they are fully automatic and  do not require any
          tuning parameters (as the shrinkage intensity is analytically estimated from the data), and
    \item they assume nothing about the underlying distributions, except for the existence of
          the first two moments.}
  These properties also carry over to derived quantities, such as partial variances and
  partial correlations (\code{\link{pvar.shrink}} and  \code{\link{pcor.shrink}}).


  
  As an extra benefit, the shrinkage estimators have a form that can be \emph{very} efficiently inverted, 
  especially if the number of variables is large and the sample size is small.   Thus, instead of 
   inverting the matrix output by \code{cov.shrink} and \code{cor.shrink} please use the functions
  \code{\link{invcov.shrink}} and \code{\link{invcor.shrink}}, respectively.
  
}
\value{
  \code{var.shrink} returns a vector with estimated variances.

  \code{cov.shrink} returns a covariance matrix.
  
  \code{cor.shrink} returns the corresponding correlation matrix. 
}
\author{
  Juliane Sch\"afer,
  Rainer Opgen-Rhein,
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
\seealso{\code{\link{invcov.shrink}}, \code{\link{pcor.shrink}}, \code{\link{cor2pcor}}}

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


# estimate covariance matrix
s1 = cov(X)
s2 = cov.shrink(X)


# squared error
sum((s1-sigma)^2)
sum((s2-sigma)^2)


# compare positive definiteness
is.positive.definite(sigma)
is.positive.definite(s1)
is.positive.definite(s2)


# compare ranks and condition
rank.condition(sigma)
rank.condition(s1)
rank.condition(s2)

# compare eigenvalues
e0 = eigen(sigma, symmetric=TRUE)$values
e1 = eigen(s1, symmetric=TRUE)$values
e2 = eigen(s2, symmetric=TRUE)$values
m = max(e0, e1, e2)
yl = c(0, m)

par(mfrow=c(1,3))
plot(e1,  main="empirical")
plot(e2,  ylim=yl, main="full shrinkage")
plot(e0,  ylim=yl, main="true")
par(mfrow=c(1,1))

}
\keyword{multivariate}
