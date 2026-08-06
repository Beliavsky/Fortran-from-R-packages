\name{YieldCurve-package}
\alias{YieldCurve-package}
\alias{YieldCurve}
\docType{package}
\title{
Modelling and estimation of the yield curve
}
\description{
Modelling the yield curve with some parametric models. 
The models implemented are: Nelson-Siegel, Diebold-Li and Svensson. 
The package also includes the data of the term structure of interest rate of Federal Reserve Bank
and European Central Bank.
}
\details{
\tabular{ll}{
Package: \tab YieldCurve\cr
Type: \tab Package\cr
Version: \tab 5\cr
Date: \tab 2022-09-30\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
DieboldLi
}
\author{
Sergio Salvino Guirreri

Maintainer: Sergio Salvino Guirreri <sergioguirreri@gmail.com>

}
\references{
Diebold, F.X. and Li, C. (2006), Forecasting the Term Structure of Government Bond Yields, \emph{Journal of Econometrics}, \bold{130}, 337-364. 

Diebold, F.X., Ji, L. and Li, C. (2006), A Three-Factor Yield Curve Model: Non-Affine Structure, Systematic Risk Sources, and Generalized Duration, in L.R. Klein (ed.), \emph{Long-Run Growth and Short-Run Stabilization: Essays in Memory of Albert Ando}. Cheltenham, U.K.: Edward Elgar, 240-274.

Nelson, C.R., and A.F. Siegel (1987), Parsimonious Modeling of Yield Curve, \emph{The Journal of Business}, \bold{60}, 473-489.

Svensson, L.E. (1994), Estimating and Interpreting Forward Interest Rates: Sweden 1992-1994, \emph{IMF Working Paper}, \bold{WP/94/114}.

}

\keyword{models}
\keyword{htest}


\examples{
### Nelson.Siegel function and Fed data-set ###
data(FedYieldCurve)
rate.Fed = first(FedYieldCurve,'5 month')
maturity.Fed <- c(3/12, 0.5, 1,2,3,5,7,10)
NSParameters <- Nelson.Siegel( rate= rate.Fed, maturity=maturity.Fed )
y <- NSrates(NSParameters[5,], maturity.Fed)
plot(maturity.Fed,rate.Fed[5,],main="Fitting Nelson-Siegel yield curve", type="o")
lines(maturity.Fed,y, col=2)
legend("topleft",legend=c("observed yield curve","fitted yield curve"),
col=c(1,2),lty=1)

### Svensson function and ECB data-set ###
data(ECBYieldCurve)
rate.ECB = ECBYieldCurve[1:5,]
maturity.ECB = c(0.25,0.5,seq(1,30,by=1))
SvenssonParameters <- Svensson(rate.ECB, maturity.ECB)
Svensson.rate <- Srates( SvenssonParameters ,maturity.ECB,"Spot")

plot(maturity.ECB, rate.ECB[5,],main="Fitting Svensson yield curve", type="o")
lines(maturity.ECB, Svensson.rate[5,], col=2)
legend("topleft",legend=c("observed yield curve","fitted yield curve"),
col=c(1,2),lty=1)
}

