\name{CG.test}
\alias{CG.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Testing survival difference of two groups via the CG estimators
}
\description{
Testing survival difference of two prognostic groups separated by a prognostic index (PI).
 Survival probabilities are computed by the CG estimators (Yeh, et al. 2023).
}
\usage{
CG.test(t.vec,d.vec,PI,cutoff=median(PI),alpha=2,
copula=CG.Clayton,S.plot=TRUE,N=10000,mark.time=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t.vec}{Vector of survival times (time to either death or censoring)}
  \item{d.vec}{Vector of censoring indicators, 1=death, 0=censoring}
  \item{PI}{Vector of real numbers (the values of a prognostic index)}
  \item{cutoff}{A number determining the cut-off value of a prognostic index}
  \item{alpha}{Copula parameter}
  \item{copula}{Copula function: "CG.Clayton","CG.Gumbel" or "CG.Frank"}
  \item{S.plot}{If TRUE, the survival curve is displayed}
  \item{N}{The number of permutations}
  \item{mark.time}{If TRUE, then curves are marked at each censoring time}
}
\details{
Two-sample comparison based on estimated survival functions
via copula-graphic estimators under dependent censoring.
The D statistic (the mean vertical difference betewen two
estimated survival functions) is used for testing the null
hypothesis of no difference in survival.
See Yeh et al.(2023) for details.
}
\value{
\item{test}{Testing the difference of two survival functions}
\item{Good}{Good prognostic group defined by PI<=c}
\item{Poor}{Poor prognostic group defined by PI>c}
}
\references{
Emura T, Chen YH (2018). Analysis of Survival Data with Dependent Censoring,
Copula-Based Approaches, JSS Research Series in Statistics, Springer, Singapore.

Rivest LP, Wells MT (2001). A Martingale Approach to the Copula-graphic Estimator for the
Survival Function under Dependent Censoring, J Multivar Anal; 79: 138-55.

Yeh CT, Liao GY, Emura T (2023). Sensitivity analysis for survival prognostic prediction
with gene selection: a copula method for dependent censoring, Biomedicines 11(3):797.
}
\author{Takeshi Emura, Pauline Baur}

\examples{
t.vec=c(1,3,5,4,7,8,10,13)
d.vec=c(1,0,0,1,1,0,1,0)
PI=c(8,7,6,5,4,3,2,1)

CG.test(t.vec,d.vec,PI,copula=CG.Clayton,alpha=18,N=100)
CG.test(t.vec,d.vec,PI,copula=CG.Gumbel,alpha=2,N=100)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dependent censoring }
\keyword{ copula }
\keyword{ copula-graphic estimator }
