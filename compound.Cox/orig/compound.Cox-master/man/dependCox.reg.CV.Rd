\name{dependCox.reg.CV}
\alias{dependCox.reg.CV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cox regression under dependent censoring.
}
\description{
This function performs estimation and significance testing for survival data under a copula-based dependent censoring model proposed in
Emura and Chen (2016). The dependency between the failure and censoring times is modeled via the Clayton copula.
The method is based on the semiparametric maximum likelihood estimation, where the association parameter is estimated by maximizing
the cross-validated c-index (see Emura and Chen 2016 for details).
}
\usage{
dependCox.reg.CV(t.vec, d.vec, X.mat, K = 5, G = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t.vec}{A vector of survival times (time-to-death or censoring)}
  \item{d.vec}{A vector of censoring indicators, 1=death, 0=censoring}
  \item{X.mat}{An (n*p) matrix of covariates, where n is the sample size and p is the number of covariates}
  \item{K}{The number of cross-validation folds}
  \item{G}{The number of grids to optimize c-index (c-index is computed for G different values of copula parameters)}
}
\details{
The Clayton model yields positive association between survival time and censoring time with Kendall's tau being equal to 
alpha/(alpha+2), where alpha > 0 is a copula parameter. The independence copula corresponds to alpha = 0. 

If the number of covariates p is large (p>=100), the computational time becomes very long. 
We suggest using "uni.selection" to reduce the number such that p<100. 

If the number of grids G is large, the computational time becomes very long. Please take 5<=G<=20.
}
\value{
\item{beta }{The estimated regression coefficients}
\item{SE }{The standard errors for the estimated regression coefficients}
\item{Z }{The Z-values for testing the null hypothesis of "beta=0" (the Wald test)}
\item{P }{The P-values for testing the null hypothesis of "beta=0" (the Wald test)}
\item{alpha }{The estimated copula parameter by optimizing c-index}
\item{c_index }{The optimized value of c_index}
}
\references{
Emura T, Chen YH (2016). Gene Selection for Survival Data Under Dependent Censoring: 
a Copula-based Approach, Stat Methods Med Res 25(No.6): 2840-57
}
\author{Takeshi Emura}


\examples{
### Reproduce Section 5 of Emura and Chen (2016) ###
data(Lung)
temp=Lung[,"train"]==TRUE
t.vec=Lung[temp,"t.vec"]
d.vec=Lung[temp,"d.vec"]
X.mat=as.matrix(Lung[temp,-c(1,2,3)])
#dependCox.reg.CV(t.vec,d.vec,X.mat,G=20) # time-consuming process #
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dependent censoring }
\keyword{ copula }
\keyword{ univariate Cox regression }% __ONLY ONE__ keyword per line
