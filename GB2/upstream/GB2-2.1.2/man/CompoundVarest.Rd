\name{CompoundVarest}
\Rdversion{2.1}
\alias{CompoundVarest}
\alias{scoreU.cgb2}
\alias{varscore.cgb2}
\alias{desvar.cgb2}
\alias{hess.cgb2}
\alias{vepar.cgb2}
\alias{derivind.cgb2}
\alias{veind.cgb2}
\title{
Variance Estimation of the Compound GB2 Distribution
}
\description{
Calculation of variance estimates of the parameters of the compound GB2 distribution and of the estimated compound GB2 indicators under cluster sampling.
}
\usage{
scoreU.cgb2(fac, pl)
varscore.cgb2(U, w=rep(1,dim(U)[1]))
desvar.cgb2(data=data, U=U, ids=NULL, probs=NULL, strata = NULL, variables = NULL, 
fpc=NULL, nest = FALSE, check.strata = !nest, weights=NULL, pps=FALSE, 
variance=c("HT","YG"))
hess.cgb2(U, pl, w=rep(1,dim(U)[1]))
vepar.cgb2(ml, Vsc, hess)
derivind.cgb2(shape1, scale, shape2, shape3, pl0, pl, prop=0.6, decomp="r")
veind.cgb2(Vpar, shape1, scale, shape2, shape3, pl0, pl, decomp="r")
}
\arguments{
  \item{fac}{numeric; matrix of Gamma factors (output of \code{fac.cgb2}.}
  \item{pl}{numeric; a vector of fitted mixture probabilities. Sums to one. If \code{pl} is equal to \code{pl0}, we obtain the GB2 distribution.}
  \item{U}{numeric; vector of scores. Output of the \code{scoreU.cgb2} function.}
  \item{w}{numeric; vector of some extrapolation weights. By default \code{w} is a vector of 1.}
  \item{data}{dataset containing the design information per unit.}
  \item{ids, probs, strata, variables, fpc, nest, check.strata, weights, pps, variance}{parameters of \code{\link[survey]{svydesign}}.}
  \item{ml}{numeric; output of the \code{ml.cgb2} function. A list with two components. First component: estimated mixture probabilities. Second component: list containing the output of \code{optim}.}
  \item{Vsc}{numeric; 4 by 4 matrix.}
  \item{hess}{numeric; Hessian (bread) for the sandwich variance estimate.}
  \item{shape1, scale ,shape2, shape3}{numeric; positive parameters of the GB2 distribution.}
  \item{pl0}{numeric; a vector of initial proportions defining the number of components and the weight of each component density in the decomposition. Sums to one.}
  \item{prop}{numeric; proportion (in general is set to 0.6).}
  \item{decomp}{string; specifying if the decomposition of the GB2 is done with respect to the right tail ("r") or the left tail ("l") of the distribution.}
  \item{Vpar}{numeric; 4 by 4 matrix. Output of the function \code{vepar.cgb2}.}
}
\details{
Function \code{scoreU.cgb2} calculates the \eqn{N \times (L-1)} matrix of scores \code{U} is defined as
\deqn{U(k,\ell)=p_{\ell} \left( \frac{F(k,\ell)}{\sum_{j=1}^L p_{j}\, F(k,j)} - 1\right),} where \eqn{p_\ell, \ell=1,..,L} is the vector of fitted mixture probabilities and \eqn{F} is the \eqn{N \times L}
matrix of gamma factors, output of \code{fg.cgb2}.
The linearized scores are the columns of \code{U}. They serve to compute the linearization approximation of the covariance matrix of the parameters \eqn{v_\ell=\log(p_\ell/p_L), 
\ell=1,...,L-1}. Function \code{varscore.cgb2} calculates the middle term of the sandwich variance estimator, that is the \eqn{((L-1) \times (L-1))} 
estimated variance-covariance matrix of the \eqn{(L-1)} weighted sums of  the columns of \code{U},  without design information. 
\code{desvar.cgb2} calculates the design-based variance-covariance matrix of the \eqn{(L-1)} weighted sums of the columns of \code{U}, invoking \code{\link[survey]{svydesign}}
and \code{\link[survey]{svytotal}} of package \pkg{survey}. \code{hess.cgb2} calculates the Hessian (\eqn{(L-1) \times (L-1))} matrix of second derivatives of the pseudo-log-likelihood with 
respect to the parameters \eqn{v_\ell}). It should be negative definite. If not, the maximum likelihood estimates are spurious. 
\code{vepar.cgb2} calculates the sandwich covariance matrix estimate of the vector of parameters \code{v}.
\code{veind.cgb2} calculates estimates, standard error, covariance and correlation matrices of the indicators under the compound GB2.
}
\value{
\code{scoreU.cgb2} returns a \eqn{N \times (L-1)} matrix of scores \code{U}.
\code{varscore.cgb2} returns the variance-covariance estimate of the  weighted sums of scores \code{U}, given by weighted cross products.
\code{desvar.cgb2} returns a list of two elements. The first is the output of svytotal and the second is the design-based variance-covariance matrix of the weighted sums of the scores. 
\code{hess.cgb2} returns the matrix of second derivatives of the likelihood with respect to the parameters (bread for the sandwich variance estimate). 
\code{vepar.cgb2} returns a list of five elements - [["type"]] with value "parameter", [["estimate"]] estimated parameters, [["stderr"]] corresponding standard errors, 
[["Vcov"]] variance -covariance matrix and [["Vcor"]] - correlation matrix.
\code{veind.cgb2} returns a list of five elements: [["type"]] with value "indicator", 
[["estimate"]] estimated indicators under the compound GB2, [["stderr"]] corresponding standard errors, [["Vcov"]] variance -covariance matrix and [["Vcor"]] - correlation matrix.
}
\references{
Davison, A. (2003),
\emph{Statistical Models}.
Cambridge University Press.

Freedman, D. A. (2006),
On The So-Called "Huber Sandwich Estimator" and "Robust Standard Errors".
\emph{The American Statistician}, \bold{60}, 299--302.

Graf, M., Nedyalkova, D., Muennich, R., Seger, J. and Zins, S. (2011)
 AMELI Deliverable 2.1: Parametric Estimation of Income Distributions and
 Indicators of Poverty and Social Exclusion.
 \emph{Technical report}, AMELI-Project.
 
Pfeffermann, D. and Sverchkov, M. Yu. (2003),
Fitting Generalized Linear Models under Informative Sampling.
In, Skinner, C.J. and Chambers, R.L. (eds.). 
\emph{Analysis of Survey Data}, chapter 12, 175--195.
Wiley, New York.
}
\author{
Monique Graf and Desislava Nedyalkova
}
\examples{
\dontrun{
# Example (following of example in CompoundFit)

# Estimated mixture probabilities:
(pl.hat <- estim[[1]])

# scores per unit
U <- scoreU.cgb2(fac, pl.hat)  

# Conditional variances given a,b,p,q:

# 1. Variance of sum of scores:
(Vsc <- t(U) %*% U)
(Vsc <- varscore.cgb2(U))        

# 2. sandwich variance-covariance matrix estimate of (v_1,v_2):
(hess <-  hess.cgb2(U,pl.hat))
(Parameters <- vepar.cgb2(estim, Vsc, hess))

# 3. Theoretical indicators  (with mixture prob pl)
decomp <- "r"
(theoretical <- main.cgb2( 0.6,a,b,p,q,pl0, pl,decomp=decomp))

# Estimated indicators and conditional variances : takes a long time!
(Indic <- veind.cgb2(Parameters,a,b,p,q, pl0, pl.hat, decomp="r") )

}
}
\keyword{distribution}
