\name{CompoundAuxVarest}
\Rdversion{2.1}
\alias{CompoundAuxVarest}
\alias{scoreU.cavgb2}
\alias{scorez.cavgb2}
\alias{varscore.cavgb2}
\alias{desvar.cavgb2}
\alias{hess.cavgb2}
\alias{vepar.cavgb2}
\alias{veind.cavgb2}
\title{
Variance Estimation under the Compound GB2 Distribution Using Auxiliary Information
}
\description{
Calculation of variance estimates of the parameters of the compound GB2 distribution and of the estimated compound GB2 indicators under a complex survey design (see package \pkg{survey}).  
}
\usage{
scoreU.cavgb2(fac, z, lambda) 
scorez.cavgb2(U,z)
varscore.cavgb2(SC, w=rep(1,dim(SC)[1]))
desvar.cavgb2(data=data, SC=SC, ids=NULL, probs=NULL, strata = NULL, variables = NULL,
fpc=NULL, nest = FALSE, check.strata = !nest, weights=NULL, pps=FALSE, 
variance=c("HT","YG"))
hess.cavgb2(U, P, z, w=rep(1, dim(z)[1]))
vepar.cavgb2(ml, Vsc, hess)
veind.cavgb2(group, vepar, shape1, scale, shape2, shape3, pl0, P, decomp="r")
}
\arguments{
  \item{fac}{numeric; a matrix of Gamma factors.}
  \item{z}{numeric; a matrix of auxiliary variables.}
  \item{lambda}{numeric; a matrix of parameters.}
  \item{U}{numeric; a matrix of scores \eqn{U_{k,\ell}} (output of the \code{scoreU.cavgb2} function).}
  \item{SC}{numeric; scores, output of \code{scorez.cavgb2}.}
  \item{w}{numeric; vector of extrapolation weights. By default \code{w} is a vector of 1.}
  \item{data}{dataset containing the design information per unit}
  \item{ids, probs, strata, variables, fpc, nest, check.strata, weights, pps, variance}{parameters of \code{\link[survey]{svydesign}}.}
  %\item{U}{numeric; scores computed in \code{scoreU.cavgb2}.}
  \item{P}{numeric; matrix of mixture probabilities (output of \code{pkl.cavgb2}).}
  %\item{z}{numeric; a matrix of auxiliary variables.} 
  \item{ml}{numeric; estimated values of the vector of v's. Output of the \code{ml.cavgb2} function (the second element in the list).}
  \item{Vsc}{numeric; 4 by 4 matrix. Variance of the scores \code{SC}, computed in \code{varscore.cavgb2} or with the design information in \code{\link{desvar.cavgb2}}.}
  \item{hess}{numeric; Hessian (bread) for the sandwich variance estimate (output of \code{hess.cavgb2}).}
  \item{group}{numeric; a factor variable of the same length as the sample size giving the group membership in the special case when the auxiliary information defines group membership.}
  \item{vepar}{numeric; output of \code{vepar.cavgb2}.}
  \item{shape1, scale ,shape2, shape3}{numeric; positive parameters of the GB2 distribution.}
  \item{pl0}{numeric; a vector of initial proportions defining the number of components and the weight of each component density in the decomposition. Sums to one.}
  \item{decomp}{string; specifying if the decomposition of the GB2 is done with respect to the right tail ("r") or the left tail ("l") of the distribution.}
}
\details{
The \eqn{N \times L} matrix of fitted mixture probabilities \code{P}\eqn{=(p_{k,\ell})} depends on the \eqn{N \times I} matrix \code{z} of auxiliary variables. 
\code{P} has as many distinct rows as there are distinct rows in \code{z}. The \eqn{N \times L} matrix  of gamma factors \code{fac}\eqn{=F}, output of \code{fg.cgb2} depends on the vector of initial 
probabilities \eqn{p_{0,\ell}} only. The \eqn{N \times (L-1)} matrix of scores \code{U} is defined as \deqn{U(k,\ell)=p_{k,\ell} \left( \frac{F(k,\ell)}{\sum_{j=1}^L p_{k,j}\, F(k,j)} - 1\right).} 
The linearized scores are the columns of a \eqn{N \times I(L-1)} matrix \deqn{SC(k,\,I(\ell-1)+i)= U(k,\ell) \,  z(k,i).} 
Function \code{varscore.cavgb2} calculates the middle term of the sandwich variance estimator, that is the \eqn{(I(L-1) \times I(L-1))} estimated variance-covariance matrix of the \eqn{I(L-1)}
weighted sums of  the columns of \eqn{SC}, without design information. 
\code{desvar.cavgb2} calculates the design-based variance-covariance matrix of the \eqn{I(L-1)} weighted sums of the columns of \eqn{SC}, invoking \code{\link[survey]{svydesign}} and \code{\link[survey]{svytotal}} of package \pkg{survey}. 
\code{hess.cavgb2} calculates the Hessian (\eqn{I(L-1) \times I(L-1)} matrix of second derivatives of the pseudo-log-likelihood with respect to the parameters). It should be negative definite. 
If not, the maximum likelihood estimates are spurious. \code{vepar.cavgb2} calculates the sandwich variance estimate of the vectorized matrix of parameters \code{lambda}.
\code{veind.cavgb2} calculates estimates, std error, covariance and correlation matrices of the indicators under the compound GB2 with auxiliary variables in the particular 
case where the unique combinations of the auxiliary variables define a small number of groups. Group membership is specified by the vector \code{group} of length \eqn{N}.
}
\value{
\code{scoreU.cavgb2} returns a \eqn{N \times (L-1)} matrix of scores \code{U}. 
\code{scorez.cavgb2} returns a \eqn{N \times I(L-1)} matrix whose columns are the linearized scores \code{SC}. 
\code{varscore.cavgb2} returns the variance-covariance estimate of the weighted sums of scores \code{SC}, given by weighted cross products. 
\code{desvar.cavgb2} returns a list of two elements. The first is the output of svytotal and the second is the design-based variance-covariance matrix of the 
weighted sums of the scores SC. 
\code{hess.cavgb2} returns the matrix of second derivatives of the likelihood with respect to the parameters (bread for the sandwich variance estimate). 
\code{vepar.cavgb2} returns a list of five elements - [["type"]] with value "parameter", [["estimate"]] estimated parameters, [["stderr"]] corresponding standard errors, [["Vcov"]] variance -covariance matrix and [["Vcor"]] - correlation matrix. 
\code{veind.cavgb2} returns a list of five elements: [["type"]] with value "indicator", followed by a list with as many arguments as \code{length(levels(group))}. Each argument is itself a list with 5
arguments: [["group"]] group name, [["estimate"]] estimated indicators under the compound GB2, [["stderr"]] corresponding standard errors, [["Vcov"]] variance -covariance matrix and [["Vcor"]] - correlation matrix.
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
# Example (following of example in CompoundAuxFit)

# Scores U 
U <- scoreU.cavgb2(facgl, z, lambdafitl)

# Scores multiplied by z
SC <- scorez.cavgb2(U,z)

# Naive variance estimate of sum of scores
(Vsc <- varscore.cavgb2(SC,w=pwa))

# Design based variance of sum of scores
(desv <- desvar.cavgb2(data=ns,SC=SC,id=~hid,strata=~region,weights=~pwa))

# Hessian
hess <- hess.cavgb2(U,pglfitl,z,w=pwa) 

# 1. Sandwich variance-covariance matrix estimate of parameters using Vsc:
Param1 <- vepar.cavgb2(fitcml,Vsc, hess)
Param1

# 2. Sandwich variance-covariance matrix estimate of parameters using 
# the design variance:
Param2 <- vepar.cavgb2(fitcml,desv$Vtheta, hess)
Param2 

# 3. Indicators and conditional variances : takes a long time!
(Indic <- veind.cavgb2(group,Param2 ,agl.fit,bgl.fit,pgl.fit,qgl.fit, 
                       pl0, pglfitl, decomp="l") )
}
}
\keyword{distribution}
