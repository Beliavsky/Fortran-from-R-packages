\name{CG.Clayton}
\alias{CG.Clayton}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Copula-graphic estimator under the Clayton copula.
}
\description{
This function computes the copula-graphic (CG) estimator (Rivest & Wells 2001) of a survival function under the Clayton copula.
}
\usage{
CG.Clayton(t.vec, d.vec, alpha, S.plot = TRUE, S.col = "black")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t.vec}{Vector of survival times (time to either death or censoring)}
  \item{d.vec}{Vector of censoring indicators, 1=death, 0=censoring}
  \item{alpha}{Association parameter that is related to Kendall's tau through "tau= alpha/(alpha+2)"}
  \item{S.plot}{If TRUE, the survival function is displayed}
  \item{S.col}{Color of the survival function in the plot}
}
\details{
The CG estimator is a variant of the Kaplan-Meier estimator for a survival function.
The CG estimator relaxes the independent censoring assumption of the KM estimator through a copula-based dependent censoring model.
The computational formula of the CG estimator is given in Appendix D of Emura et al. (2019) or Section 3.2 of Yeh et al.(2023).
The output shows the survival probabilities at given time points of "t.vec".
The input requires to specify an association parameter "alpha" of the Clayton copula (alpha>0), where alpha=0 corresponds to the independence copula.
Emura and Chen (2016, 2018) and Yeh et al.(2023) applied the CG estimator to assess survival prognosis for lung cancer patients.
}
\value{
\item{tau }{Kendall's tau (=alpha/(alpha+2))}
\item{time }{sort(t.vec)}
\item{n.risk }{the number of patients at-risk}
\item{surv }{survival probability at "time"}
}
\references{
Emura T, Matsui S, Chen HY (2019). compound.Cox: Univariate Feature Selection and Compound Covariate for Predicting Survival,
Computer Methods and Programs in Biomedicine 168: 21-37.

Emura T, Chen YH (2016). Gene Selection for Survival Data Under Dependent Censoring:
a Copula-based Approach, Stat Methods Med Res 25(No.6): 2840-57.

Emura T, Chen YH (2018). Analysis of Survival Data with Dependent Censoring,
Copula-Based Approaches, JSS Research Series in Statistics, Springer, Singapore.

Rivest LP, Wells MT (2001). A Martingale Approach to the Copula-graphic Estimator for the
Survival Function under Dependent Censoring, J Multivar Anal; 79: 138-55.

Yeh CT, Liao GY, Emura T (2023). Sensitivity analysis for survival prognostic prediction
with gene selection: a copula method for dependent censoring, Biomedicines 11(3):797.
}
\author{Takeshi Emura}

\examples{
## Example 1 (a toy example of n=8) ##
t.vec=c(1,3,5,4,7,8,10,13)
d.vec=c(1,0,0,1,1,0,1,0)
CG.Clayton(t.vec,d.vec,alpha=18,S.col="blue")
### CG.Clayton gives identical results with the Kaplan-Meier estimator with alpha=0 ###
CG.Clayton(t.vec,d.vec,alpha=0,S.plot=FALSE)$surv
survfit(Surv(t.vec,d.vec)~1)$surv

## Example 2 (Analysis of the lung cancer data) ##
data(Lung) # read the data
t.vec=Lung[,"t.vec"]
d.vec=Lung[,"d.vec"]
x.vec=Lung[,"MMP16"] # the gene associated with survival (Emura and Chen 2016, 2018) #
Poor=x.vec>median(x.vec) ## Indicator of poor survival
Good=x.vec<=median(x.vec) ## Indicator of good survival

par(mfrow=c(1,2))
###### Predicted survival curves via the CG estimator #####
t.good=t.vec[Good]
d.good=d.vec[Good]
CG.Clayton(t.good,d.good,alpha=18,S.plot=TRUE,S.col="blue")

t.poor=t.vec[Poor]
d.poor=d.vec[Poor]
CG.Clayton(t.poor,d.poor,alpha=18,S.plot=TRUE,S.col="red")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ dependent censoring }
\keyword{ copula }
\keyword{ copula-graphic estimator }
