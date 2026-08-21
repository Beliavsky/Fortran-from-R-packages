\name{CompoundFit}
\Rdversion{2.1}
\alias{CompoundFit}
\alias{vofp.cgb2}
\alias{pofv.cgb2}
\alias{logl.cgb2}
\alias{scores.cgb2}
\alias{ml.cgb2}
\title{
Fitting the Compound Distribution based on the GB2 by the Method of Maximum Likelihood Estimation
}
\description{
Calculates the log-likelihood, the score functions of the log-likelihood, the weighted mean of scores, and fits the parameters of the Compound Distribution based on the GB2.
}
\usage{
vofp.cgb2(pl)
pofv.cgb2(vl)
logl.cgb2(fac, pl, w=rep(1, dim(fac)[1]))
scores.cgb2(fac, pl, w=rep(1, dim(fac)[1]))
ml.cgb2(fac, pl0, w=rep(1, dim(fac)[1]), maxiter=100, fnscale=length(w)) 
}
\arguments{
  \item{pl0}{numeric; vector of initial proportions defining the number of components and the weight of each component density in the decomposition. Sums to one.}
  \item{pl}{numeric; vector of fitted proportions. Sums to one. If \code{pl} is equal to \code{pl0}, we obtain the GB2 distribution.}
  \item{fac}{numeric; matrix of Gamma factors (output of \code{fac.cgb2}.}
  \item{vl}{numeric; vector of parameters. Its length is equal to the length of \code{pl} - 1.}
  \item{w}{numeric; vector of weights of length the number of rows of the matrix \code{fac}. By default \code{w} is a vector of 1.}
  \item{maxiter}{numeric; maximum number of iterations to perform. By default \code{maxiter} = 100.}
  \item{fnscale}{numeric; an overall scaling parameter used in the function \code{\link[stats]{optim}}. By default it is equal to the length of the vector of weights \code{w}.}
}
\details{
There are only \eqn{L-1} parameters to estimate, because the probabilities \eqn{p_\ell} sum to 1 (L is the dimension of the vector of probailities \eqn{p_\ell}). Knowing this, we change the parameters \eqn{p_\ell} to \eqn{v_\ell=log(p_\ell/p_L), \ \ell= 1, ..., L-1}. This calculation is done through the function \code{vofp.cgb2}. \code{pofv.cgb2} calculates the \eqn{p_\ell} in function of the given \eqn{v_\ell}.
We express the log-likelihood as a weighted mean of \eqn{log(f) = log(\sum(p_\ell f_\ell)}, evaluated at the data points, where \eqn{f} is the GB2 compound density. If the weights are not available, then we suppose that \code{w} \eqn{= 1}. Analogically, the scores are obtained as weighted sums of the first derivatives of the log-likelihood, with respect to the parameters \eqn{v_\ell, \ \ell=1, ..., L-1}, evaluated at the data points.  Function \code{ml.cgb2} performs maximum likelihood estimation through the general-purpose optimization function \code{optim} from package \code{stats}. The considered method of optimization is BFGS.
}
\value{
\code{vofp.cgb2} returns a vector of length \eqn{L-1}, where \eqn{L} is the length of the vector \eqn{p_\ell}. 
\code{pofv.cgb2} returns a vector of length \eqn{\ell}. 
\code{logl.cgb2} returns the value of the pseudo log-likelihood. 
\code{scores.cgb2} returns a vector of the weighted mean of the scores of length \eqn{L-1}. 
\code{ml.cgb2} returns a list containing two objects - the vector of fitted proportions \eqn{\hat{p_\ell}} and the output of the BFGS fit.
}
\seealso{
\code{\link[stats]{optim}}
}
\author{
Monique Graf and Desislava Nedyalkova
}
\examples{
\dontrun{
# GB2 parameters:
a <- 4
b <- 1950
p <- 0.8
q <- 0.6

# Proportions defining the component densities:
pl0 <- rep(1/3,3)

# Mixture probabilities
pl <- c(0.1,0.8,0.1)

# Random generation:
n <- 10000
set.seed(12345)
x <- rcgb2(n,a,b,p,q,pl0,pl,decomp="l")

# Factors in component densities
fac <- fg.cgb2(x,a,b,p,q, pl0,decomp="l")

# Estimate the mixture probabilities:
estim <- ml.cgb2(fac,pl0)

# estimated mixture probabilities:
estim[[1]]
#[1] 0.09724319 0.78415797 0.11859883
}
}
\keyword{distribution}


