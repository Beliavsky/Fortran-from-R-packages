\name{cindex.CV}
\alias{cindex.CV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross-validated c-index for measuring the 
predictive accuracy of a prognostic index under a copula-based dependent censoring model.
}
\description{
This function calculates the cross-validated c-index (concordance index) for measuring the 
predictive accuracy of a prognostic index under a copula-based dependent censoring model. 
Here the prognostic index is calculated as a compound covariate predictor based on the 
univariate Cox regression estimates. The expression and details are given in Section 3.2 of Emura and Chen (2016). 
The association between survival time and censoring time is modeled via the Clayton copula.
}
\usage{
cindex.CV(t.vec, d.vec, X.mat, alpha, K = 5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t.vec}{Vector of survival times (time to death or time to censoring, whichever comes first)}
  \item{d.vec}{Vector of censoring indicators, 1=death, 0=censoring}
  \item{X.mat}{n by p matrix of covariates, where n is the sample size and p is the number of covariates}
  \item{alpha}{Association parameter of the Clayton copula; Kendall's tau = alpha/(alpha+2)}
  \item{K}{The number of cross-validation folds (K=5 is the defailt)}
}
\details{
Currently, only the Clayton copula is implemented for modeling association between survival time 
and censoring time. The Clayton model yields positive association between survival time and censoring 
time with the Kendall's tau being equal to alpha/(alpha+2), 
where alpha > 0. The independent copula corresponds to alpha = 0.

If the number of covariates p is large (e.g., p>=100), the computational time becomes very long. 
Pre-filtering for covariates is recommended to reduce p. 
}
\value{
\item{concordant }{Cross-validated c-index}
}

\references{
Emura T, Chen YH (2016). Gene Selection for Survival Data Under Dependent Censoring: 
a Copula-based Approach, Stat Methods Med Res 25(No.6): 2840-57.
}
\author{Takeshi Emura}

\examples{

n=25 ### sample size ###
p=3  ### the number of covariates ###
set.seed(1)
T=rexp(n) ### survival time
U=rexp(n) ### censoring time
t.vec=pmin(T,U) ### minimum of survival time and censoring time
d.vec=as.numeric( c(T<=U) ) ### censoring indicator
X.mat=matrix(runif(n*p),n,p) ### covariates matrix

cindex.CV(t.vec,d.vec,X.mat,alpha=2) ### alpha=2 corresponds to Kendall's tau=0.5

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cross-validation }
\keyword{ copula }
\keyword{ c-index }
