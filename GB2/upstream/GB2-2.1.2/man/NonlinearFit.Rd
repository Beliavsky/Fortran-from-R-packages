\name{NonlinearFit}
\Rdversion{2.1}
\alias{NonlinearFit}
\alias{nlsfit.gb2}
\title{
Fitting the GB2 by Minimizing the Distance Between a Set of Empirical Indicators and Their GB2 Expressions
}
\description{
Fitting the parameters of the GB2 distribution by optimizing the squared weighted distance between a set of empirical indicators, i.e. the median, the ARPR, the RMPG, the QSR and the Gini coefficient,
and the GB2 indicators using nonlinear least squares (function \code{nls} from package \code{stats}).
}
\usage{
nlsfit.gb2(med, ei4, par0=c(1/ei4[4],med,1,1), cva=1, bound1=par0[1]*max(0.2,1-2*cva), 
bound2=par0[1]*min(2,1+2*cva), ei4w=1/ei4)
}
\arguments{
  \item{med}{numeric; the empirical median.}
  \item{ei4}{numeric; the values of the empirical indicators.}
  \item{par0}{numeric; vector of initial values for the GB2 parameters \eqn{a, b, p} and \eqn{q}. The default is to take \eqn{a} equal to the inverse of the empirical Gini coefficient, \eqn{b} equal to the empirical median and \eqn{p = q = 1}.}
  \item{cva}{numeric; the coefficient of variation of the ML estimate of the parameter \eqn{a}. The default value is 1.}
  \item{bound1, bound2}{numeric; the lower and upper bounds for the parameter \eqn{a} in the algorithm. The default values are \eqn{0.2*a_{0}} and \eqn{2*a_{0}}, where \eqn{a_{0}} is the initial value of the parameter \eqn{a}.}
  \item{ei4w}{numeric; vector of weights of to be passed to the \code{nls} function. The default values are the inverse of the empirical indicators.}
}
\details{
 We consider the following set of indicators \eqn{A = (median, ARPR, RMPG, QSR, Gini)} and their corresponding GB2 expressions \eqn{A_{GB2}}. We fit the parameters of the GB2 in two consecutive steps. In the first step, we use the set of indicators, excluding the median, and their corresponding expressions in function of \eqn{a}, \eqn{ap} and \eqn{aq}. The bounds for \eqn{a} are defined in function of the coefficient of variation of the fitted parameter \eqn{\hat(a)}. The nonlinear model that is passed to \code{nls} is given by:
  \deqn{\sum_{i=1}^4 c_i (A_{empir,i}-A_{GB2,i}(a,ap,aq))^2,} where the weights \eqn{c_i} take the differing scales into account and are given by the vector \code{ei4w}. 
 \eqn{ap} and \eqn{aq} are bounded so that the constraints for the existence of the moments of the GB2 distribution and the excess for calculating the Gini coefficient are fulfilled, 
 i.e. \eqn{ap \ge 1} and \eqn{aq \ge 2}.  In the second step, only the the parameter \eqn{b} is estimated, optimizing the weighted square difference between the empirical median and the GB2 median in function of the already obtained NLS parameters \eqn{a, p} and \eqn{q}. 
}
\value{
\code{nlsfit.gb2} returns a list of three values: the fitted GB2 parameters, the first fitted object and the second fitted object.
}
\seealso{
\code{\link[stats]{nls}}, \code{\link{Thomae}}, \code{\link{moment.gb2}}
}
\author{
Monique Graf and Desislava Nedyalkova
}
\examples{
# Takes long time to run, as it makes a call to the function ml.gb2

\dontrun{
library(laeken)
data(eusilc)

# Personal income
inc <- as.vector(eusilc$eqIncome)

# Sampling weights
w <- eusilc$rb050

# Data set
d <- data.frame(inc, w)
d <- d[!is.na(d$inc),]
   
# Truncate at 0
d <- d[d$inc > 0,]
inc <- d$inc
w   <- d$w

# ML fit, full log-likelihood
fitf <- ml.gb2(inc, w)$opt1

# Estimated parameters
af <- fitf$par[1]
bf <- fitf$par[2]
pf <- fitf$par[3]
qf <- fitf$par[4]

gb2.par <- c(af, bf, pf, qf)

# Empirical indicators
indicEMP  <- main.emp(inc, w)
indicEMP <- c(indicEMP[1],indicEMP[3:6])
indicE <- round(indicEMP, digits=3)

# Nonlinear fit
nn <- nlsfit.gb2(indicEMP[1,3:6],indicEMP[3:6])
an <- nn[[1]][1]
bn <- nn[[1]][2]
pn <- nn[[1]][3]
qn <- nn[[1]][4]

# GB2 indicators
indicNLS <- c(main.gb2(0.6, an, bn, pn, qn)[1], main.gb2(0.6, an, bn, pn, qn)[3:6])
indicML <- c(main.gb2(0.6, af, bf, pf, qf)[1], main.gb2(0.6, af, bf, pf, qf)[3:6])
indicN <- round(indicNLS, digits=3)
indicM <- round(indicML, digits=3)

# Likelihoods
nlik <- loglp.gb2(inc, an, bn, pn, qn, w)
mlik <- loglp.gb2(inc, af, bf, pf, qf, w)

# Results
type=c("Emp. est", "NLS", "ML full")
results <- data.frame(type=type,
        median=c(indicE[1], indicN[1], indicM[1]),
        ARPR=c(indicE[2], indicN[2], indicM[2]),
        RMPG=c(indicE[3], indicN[3], indicM[3]),
        QSR =c(indicE[4], indicN[4], indicM[4]),
        GINI=c(indicE[5], indicN[5], indicM[5]),
        likelihood=c(NA, nlik, mlik),
        a=c(NA, an, af), b=c(NA, bn, bf) ,p=c(NA, pn, pf), q=c(NA, qn, qf))

}
}
\keyword{
distribution
}

