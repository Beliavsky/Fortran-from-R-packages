\name{invcov.shrink}
\alias{invcov.shrink}
\alias{invcor.shrink}


\title{Fast Computation of the Inverse of the Covariance and Correlation Matrix}
\description{
  The functions \code{invcov.shrink} and \code{invcor.shrink} implement an
  algorithm to \emph{efficiently} compute 
  the inverses of shrinkage estimates of covariance (\code{\link{cov.shrink}}) 
  and correlation (\code{\link{cor.shrink}}).
}
\usage{
invcov.shrink(x, lambda, lambda.var, w, verbose=TRUE)
invcor.shrink(x, lambda, w, verbose=TRUE)
}
\arguments{
  \item{x}{a data matrix}
  \item{lambda}{the correlation shrinkage intensity (range 0-1).
                If \code{lambda} is not specified (the default) it is estimated
		using an analytic formula from Sch\"afer and Strimmer (2005) 
		- see \code{\link{cor.shrink}}.  
		For \code{lambda=0} the empirical correlations are recovered.}
  \item{lambda.var}{the variance shrinkage intensity (range 0-1). 
                If \code{lambda.var} is not specified (the default) it is estimated
	        using an analytic formula from Sch\"afer and Strimmer (2005) 
                - see \code{\link{var.shrink}}.  
	        For \code{lambda.var=0} the empirical variances are recovered.}	 
  \item{w}{optional: weights for each data point - if not specified uniform weights are assumed
            (\code{w = rep(1/n, n)} with \code{n = nrow(x)}).}	
  \item{verbose}{output status while computing (default: TRUE)}
}
\details{
Both \code{invcov.shrink} and \code{invcor.shrink} rely on 
\code{\link{powcor.shrink}}.  This allows to compute the inverses in
a very efficient fashion (much more efficient than directly inverting
the matrices - see the example).
}
\value{
  \code{invcov.shrink} returns the inverse of the output from \code{\link{cov.shrink}}. 
   
  \code{invcor.shrink} returns the inverse of the output from \code{\link{cor.shrink}}. 
}
\author{
  Juliane Sch\"afer  
  and Korbinian Strimmer (\url{https://strimmerlab.github.io}).
  
}
\references{
  Sch\"afer, J., and K. Strimmer. 2005.  A shrinkage approach to large-scale
  covariance estimation and implications for functional genomics. 
  Statist. Appl. Genet. Mol. Biol. \bold{4}:32.
  <DOI:10.2202/1544-6115.1175>
}
\seealso{\code{\link{powcor.shrink}}, \code{\link{cov.shrink}}, \code{\link{pcor.shrink}}, \code{\link{cor2pcor}}}

\examples{
# load corpcor library
library("corpcor")

# generate data matrix
p = 500
n = 10
X = matrix(rnorm(n*p), nrow = n, ncol = p)

lambda = 0.23  # some arbitrary lambda

# slow
system.time(
  (W1 =  solve(cov.shrink(X, lambda)))
)

# very fast
system.time(
  (W2 = invcov.shrink(X, lambda))
)

# no difference
sum((W1-W2)^2)
}
\keyword{multivariate}
