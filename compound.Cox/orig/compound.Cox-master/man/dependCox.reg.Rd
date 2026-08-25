\name{dependCox.reg}
\alias{dependCox.reg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Univariate Cox regression under dependent censoring.
}
\description{
This function performs univariate Cox regression under dependent censoring,
 where dependence between survival time and censoring time is modeled via the Clayton copula (Emura and Chen 2016).
}
\usage{
dependCox.reg(t.vec, d.vec, X.vec, alpha, var = TRUE, censor.reg=FALSE, baseline=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t.vec}{A vector of survival times (time-to-death or censoring)}
  \item{d.vec}{A vector of censoring indicators, 1=death, 0=censoring}
  \item{X.vec}{A vector of covariates (multiple covariates are not allowed)}
  \item{alpha}{An copula parameter (Kendall's tau = alpha/(alpha+2)}
  \item{var}{If TRUE, the standard deviations are given (use FALSE to reduce the computational cost)}
  \item{censor.reg}{If TRUE, show the fitted results for both survival and censoring models}
  \item{baseline}{If TRUE, show the cumulative baseline hazards at the values of "t.vec"}
}
\details{
The Clayton model yields positive association between survival time and censoring time with Kendall's tau being equal to 
alpha/(alpha+2), where alpha > 0 is a copula parameter. The independence copula corresponds to alpha = 0. 
}
\value{
\item{beta }{The estimated regression coefficient}
\item{SE }{The standard error for the estimated regression coefficient}
\item{Z }{The Z-value for testing the null hypothesis of "beta=0" (the Wald test)}
\item{P }{The P-value for testing the null hypothesis of "beta=0" (the Wald test)}
}
\references{
Emura T, Chen YH (2016). Gene Selection for Survival Data Under Dependent Censoring: 
a Copula-based Approach, Stat Methods Med Res 25(No.6): 2840-57.
}
\author{Takeshi Emura}


\examples{
### Joint Cox regression of survival and censoring ### 
data(Lung)
t.vec=Lung[,"t.vec"]# death or censoring times #
d.vec=Lung[,"d.vec"]# censoring indicators #
# 16-gene prognostic index (Emura and Chen 2016; 2018) #
X.vec=0.51*Lung[,"ZNF264"]+0.50*Lung[,"MMP16"]+
  0.50*Lung[,"HGF"]-0.49*Lung[,"HCK"]+0.47*Lung[,"NF1"]+
  0.46*Lung[,"ERBB3"]+0.57*Lung[,"NR2F6"]+0.77*Lung[,"AXL"]+
  0.51*Lung[,"CDC23"]+0.92*Lung[,"DLG2"]-0.34*Lung[,"IGF2"]+
  0.54*Lung[,"RBBP6"]+0.51*Lung[,"COX11"]+
  0.40*Lung[,"DUSP6"]-0.37*Lung[,"ENG"]-0.41*Lung[,"IHPK1"]
dependCox.reg(t.vec,d.vec,X.vec,alpha=18,censor.reg=TRUE)

temp=c(Lung[,"train"]==TRUE)
t.vec=Lung[temp,"t.vec"]
d.vec=Lung[temp,"d.vec"]
dependCox.reg(t.vec,d.vec,Lung[temp,"ZNF264"],alpha=18) 
# this reproduces Table 3 of Emura and Chen (2016) #

#### A simulation study under dependent censoring ####
beta_true=1.5 # true regression coefficient
alpha_true=2 # true copula parameter corresponding to Kendall's tau=0.5
n=150
t.vec=d.vec=X.vec=numeric(n)
set.seed(1)
for(i in 1:n){
  X.vec[i]=runif(1)
  eta=X.vec[i]*beta_true
  U=runif(1)
  V=runif(1)
  T=-1/exp(eta)*log(1-U) # Exp(eta) distribution
  W=(1-U)^(-alpha_true) # dependence produced by the Clayton copula
  C=1/alpha_true/exp(eta)*log(  1-W+W*(1-V)^(-alpha_true/(alpha_true+1))  ) # Exp(eta) distribution
  t.vec[i]=min(T,C)
  d.vec[i]=(T<=C) 
}

dependCox.reg(t.vec,d.vec,X.vec,alpha=alpha_true,var=FALSE) # faster computation by "var=FALSE"
beta_true# the above estimate is close to the true value
coxph(Surv(t.vec,d.vec)~X.vec)$coef
# this estimate is biased for the true value due to dependent censoring 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dependent censoring }
\keyword{ copula }
\keyword{ univariate Cox regression }% __ONLY ONE__ keyword per line
