\name{surv.factorial}
\alias{surv.factorial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Factorial survival analysis under dependent censoring
}
\description{
Perform factorial survival analysis under dependent censoring under an assumed copula (Emura et al. 2024).
}
\usage{
surv.factorial(t.vec,d.vec,group,copula,alpha,R=1000,t.upper=min(tapply(t.vec,group,max)),
C=NULL,S.plot=TRUE,mark.time=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t.vec}{Vector of survival times (time to either death or censoring)}
  \item{d.vec}{Vector of censoring indicators, 1=death, 0=censoring}
  \item{group}{Vector of group indicators, 1, 2, ..., d}
  \item{copula}{Copula function: "CG.Clayton","CG.Gumbel" or "CG.Frank"}
  \item{alpha}{Copula parameter}
  \item{R}{The number of Monte Carlo simulations to find the critical value of the F-test}
  \item{t.upper}{Follow-up end (default is max(t.vec))}
  \item{C}{Contrast matrix}
  \item{S.plot}{If TRUE, the survival curve is displayed}
  \item{mark.time}{If TRUE, then curves are marked at each censoring time}
}
\details{
Estimates of treatment effects and the test results are shown.
}
\value{
\item{copula.parameter}{Copula parameter}
\item{p}{Estimates of treatment effects}
\item{Var}{Variance estimates}
\item{F}{F-statistic}
\item{c.simu}{Critical value via the simulation method}
\item{c.anal}{Critical value via the analytical method}
\item{P.value}{P-value of the F-test}
}
\references{
Emura T, Ditzhaus M, Dobler D, Murotani K (2024), Factorial survival analysis for treatment effects under dependent censoring, Stat Methods Med Res 33(1):61-79.
}
\author{Takeshi Emura}

\examples{

library(survival)
data(cancer)
dat=subset(colon,etype==1)
## Treatment effects ##
t.vec=dat$time
d.vec=dat$status
trt=dat$rx
C12=matrix(c(1,-1,0),1,3,byrow=TRUE)
C13=matrix(c(1,0,-1),1,3,byrow=TRUE)
C23=matrix(c(0,1,-1),1,3,byrow=TRUE)
group=as.numeric(trt) # 1=Obs; 2=Lev; 3=Lev+5FU
# surv.factorial(t.vec,d.vec,group,alpha=2,copula=CG.Clayton)
# surv.factorial(t.vec,d.vec,group,alpha=2,copula=CG.Clayton,C=C12)
# surv.factorial(t.vec,d.vec,group,alpha=2,copula=CG.Clayton,C=C13)
# surv.factorial(t.vec,d.vec,group,alpha=2,copula=CG.Clayton,C=C23)
# alpha is a copula parameter

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dependent censoring }
\keyword{ copula }
\keyword{ copula-graphic estimator }
