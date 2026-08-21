\name{Varest}
\Rdversion{2.1}
\alias{Varest}
\alias{varscore.gb2}
\alias{vepar.gb2}
\alias{derivind.gb2}
\alias{veind.gb2}
\title{
Variance Estimation of the Parameters of the GB2 Distribution
}
\description{
Calculation of variance estimates of the estimated GB2 parameters and the estimated GB2 indicators under cluster sampling.
}
\usage{
varscore.gb2(x, shape1, scale, shape2, shape3, w=rep(1,length(x)), hs=rep(1,length(x)))
vepar.gb2(x, Vsc, shape1, scale, shape2, shape3, w=rep(1,length(x)), hs=rep(1,length(x)))
derivind.gb2(shape1, scale, shape2, shape3)
veind.gb2(Vpar, shape1, scale, shape2, shape3)
}
\arguments{
  \item{x}{numeric; vector of data values.}
  \item{Vsc}{numeric; 4 by 4 matrix.}
  \item{shape1}{numeric; positive parameter.}
  \item{scale}{numeric; positive parameter.}
  \item{shape2, shape3}{numeric; positive parameters of the Beta distribution.}
  \item{w}{numeric; vector of weights. Must have the same length as \code{x}. By default \code{w} is a vector of 1.}
  \item{hs}{numeric; vector of household sizes. Must have the same length as \code{x}. By default \code{w} is a vector of 1.}
  \item{Vpar}{numeric; 4 by 4 matrix.}
}
\details{
Knowing the first and second derivatives of \eqn{log(f)}, and using the sandwich variance estimator (see Freedman (2006)), the calculation of the variance estimates of the GB2
parameters is straightforward. \code{Vsc} is a square matrix of size the number of parameters, e.g. the estimated design variance-covariance matrix of the estimated parameters. We know that the GB2 estimates of the Laeken indicators are functions of the GB2 parameters. In this case, the variance estimates of the fitted indicators are obtained 
using the delta method. The function \code{veind.gb2} uses \code{Vpar}, the sandwich variance estimator of the vector of parameters, in order to obtain the sandwich variance estimator of the indicators. More details can be found in Graf and Nedyalkova (2011).
}
\value{
\code{varscore.gb2} calculates the middle term of the sandwich variance estimator under simple random cluster sampling. \code{vepar.gb2} returns a list of two elements: 
the estimated variance-covariance matrix of the estimated GB2 parameters and the second-order partial derivative of the pseudo log-likelihood function.
The function \code{veind.gb2} returns the estimated variance-covariance matrix of the estimated GB2 indicators.  \code{derivind.gb2} calculates the numerical derivatives of the GB2 indicators and is for internal use only. 
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
# An example of variance estimation of the GB2 parameters,
# using the dataset "eusilcP" from the R package simFrame.
# Takes long time to run

\dontrun{
library(survey)
library(simFrame)
data(eusilcP)

# Draw a sample from eusilcP
# 1-stage simple random cluster sampling of size 6000 (cluster = household)
# directly,
#s <- draw(eusilcP[, c("hid", "hsize", "eqIncome")], grouping = "hid", size = 6000)

# or setting up 250 samples, and drawing the first one.
# This sample setup can be used for running a simulation.
set.seed(12345)
scs <- setup(eusilcP, grouping = "hid", size = 6000, k = 250)
s <- draw(eusilcP[, c("region", "hid", "hsize", "eqIncome")], scs, i=1)

# The number of observations (persons) in eusilcP (58654 persons)
\dontrun{N <- dim(eusilcP)[1]} 
# The number of households in eusilcP (25000 households)
Nh <- length(unique(eusilcP$hid)) 

# Survey design corresponding to the drawn sample
sdo = svydesign(id=~hid, fpc=rep(Nh,nrow(s)), data=s)
\dontrun{summary(sdo)}

# Truncated sample (truncate at 0)
s <- s[!is.na(s$eqIncome),] 
str <- s[s$eqIncome > 0, ]
eqInc <- str$eqIncome
w <- str$.weight            

# Designs for the truncated sample
sdotr <- subset(sdo, eqIncome >0)
sddtr = svydesign(id=~hid, strata=~region, fpc=NULL, weights=~.weight, data=str)
\dontrun{summary(sdotr)}
\dontrun{summary(sddtr)}

# Fit by maximum likelihood
fit <- ml.gb2(eqInc,w)$opt1
af <- fit$par[1]
bf <- fit$par[2]
pf <- fit$par[3]
qf <- fit$par[4]
mlik <- -fit$value

# Estimated parameters and indicators, empirical indicators 
gb2.par <- round(c(af, bf, pf, qf), digits=3)
emp.ind <- main.emp(eqInc, w)
gb2.ind <- main.gb2(0.6, af, bf, pf, qf)

# Scores
scores <- matrix(nrow=length(eqInc), ncol=4)
for (i in 1:length(eqInc)){
scores[i,] <- dlogf.gb2(eqInc[i], af, bf, pf, qf)
}

# Data on households only
sh <- unique(str)
heqInc <- sh$eqIncome
hw <- sh$.weight
hhs <- sh$hsize 
hs <- as.numeric(as.vector(hhs))  

# Variance of the scores
VSC <- varscore.gb2(heqInc, af, bf, pf, qf, hw, hs)

# Variance of the scores using the explicit designs, and package survey
DV1 <- vcov(svytotal(~scores[,1]+scores[,2]+scores[,3]+scores[,4], design=sdotr))
DV2 <- vcov(svytotal(~scores[,1]+scores[,2]+scores[,3]+scores[,4], design=sddtr))

# Estimated variance-covariance matrix of the parameters af, bf, pf and qf 
VCMP <- vepar.gb2(heqInc, VSC, af, bf, pf, qf, hw, hs)[[1]]
DVCMP1 <- vepar.gb2(heqInc, DV1, af, bf, pf, qf, hw, hs)[[1]]
DVCMP2 <- vepar.gb2(heqInc, DV2, af, bf, pf, qf, hw, hs)[[1]]

\dontrun{diag(DVCMP1)/diag(VCMP)}
\dontrun{diag(DVCMP2)/diag(VCMP)}
\dontrun{diag(DV1)/diag(VSC)}
\dontrun{diag(DV2)/diag(VSC)}

# Standard errors of af, bf, pf and qf
se.par <- sqrt(diag(VCMP))
sed1.par <- sqrt(diag(DVCMP1))
sed2.par <- sqrt(diag(DVCMP2))

# Estimated variance-covariance matrix of the indicators (VCMI)
VCMI <- veind.gb2(VCMP, af, bf, pf, qf) 
DVCMI1 <- veind.gb2(DVCMP1, af, bf, pf, qf)
DVCMI2 <- veind.gb2(DVCMP2, af, bf, pf, qf)

# Standard errors and confidence intervals
varest.ind <- diag(VCMI)
se.ind <- sqrt(varest.ind)
lci.ind <- gb2.ind - 1.96*se.ind
uci.ind <- gb2.ind + 1.96*se.ind
inCI <- as.numeric(lci.ind <= emp.ind & emp.ind <= uci.ind)

# under the sampling design sdotr

varestd1.ind <- diag(DVCMI1)
sed1.ind <- sqrt(varestd1.ind)
lcid1.ind <- gb2.ind - 1.96*sed1.ind
ucid1.ind <- gb2.ind + 1.96*sed1.ind
inCId1 <- as.numeric(lcid1.ind <= emp.ind & emp.ind <= ucid1.ind)

#under the sampling design sddtr

varestd2.ind <- diag(DVCMI2)
sed2.ind <- sqrt(varestd2.ind)
lcid2.ind <- gb2.ind - 1.96*sed2.ind
ucid2.ind <- gb2.ind + 1.96*sed2.ind
inCId2 <- as.numeric(lcid2.ind <= emp.ind & emp.ind <= ucid2.ind)

#coefficients of variation .par (parameters), .ind (indicators)
cv.par <- se.par/gb2.par
names(cv.par) <- c("am","bm","pm","qm")
cvd1.par <- sed1.par/gb2.par
names(cvd1.par) <- c("am","bm","pm","qm")
cvd2.par <- sed2.par/gb2.par
names(cvd2.par) <- c("am","bm","pm","qm")

cv.ind <- se.ind/gb2.ind
cvd1.ind <- sed1.ind/gb2.ind
cvd2.ind <- sed2.ind/gb2.ind


#results 
res <- data.frame(am = af, bm = bf, pm = pf, qm = qf, lik = mlik,
  median = gb2.ind[[1]], mean = gb2.ind[[2]], ARPR = gb2.ind[[3]], 
    RMPG = gb2.ind[[4]], QSR = gb2.ind[[5]], Gini = gb2.ind[[6]],
  emedian = emp.ind[[1]], emean = emp.ind[[2]], eARPR = emp.ind[[3]], 
    eRMPG = emp.ind[[4]], eQSR = emp.ind[[5]], eGini = emp.ind[[6]],
  cva = cv.par[1], cvb = cv.par[2], cvp= cv.par[3], cvq = cv.par[4],
  cvd1a = cvd1.par[1], cvd1b = cvd1.par[2], cvd1p= cvd1.par[3], cvd1q = cvd1.par[4],
  cvd2a = cvd2.par[1], cvd2b = cvd2.par[2], cvd2p= cvd2.par[3], cvd2q = cvd2.par[4],
  cvmed = cv.ind[[1]], cvmean = cv.ind[[2]], cvARPR = cv.ind[[3]], 
  cvRMPG = cv.ind[[4]], cvQSR = cv.ind[[5]], cvGini = cv.ind[[6]],
  cvd1med = cvd1.ind[[1]], cvd1mean = cvd1.ind[[2]], cvd1ARPR = cvd1.ind[[3]], 
  cvd1RMPG = cvd1.ind[[4]], cvd1QSR = cvd1.ind[[5]], cvd1Gini = cvd1.ind[[6]],
  cvd2med = cvd2.ind[[1]], cvd2mean = cvd2.ind[[2]], cvd2ARPR = cvd2.ind[[3]], 
  cvd2RMPG = cvd2.ind[[4]], cvd2QSR = cvd2.ind[[5]], cvd2Gini = cvd2.ind[[6]])

  res <- list(parameters = data.frame(am = af, bm = bf, pm = pf, qm = qf, lik = mlik),
              cv.parameters.naive = cv.par,
              cv.parameters.design1 = cvd1.par,
              cv.parameters.design2 = cvd2.par,
		  GB2.indicators = gb2.ind,
              emp.indicators = emp.ind,
              cv.indicators.naive = cv.ind,
              cv.indicators.design1 = cvd1.ind,
              cv.indicators.design2 = cvd2.ind)
res
\dontrun{inCI}
}
}
\keyword{distribution}
