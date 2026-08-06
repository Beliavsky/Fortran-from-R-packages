\name{NSrates}
\alias{NSrates}
\title{ Interest rates of the Nelson-Siegel's model. }
\description{
Returns the interest rates by Nelson-Siegel's model.
}
\usage{
NSrates(Coeff, maturity)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Coeff}{ Vector or matrix of the beta's coefficients and lambda as the function \code{Nelson.Siegel} returns. }
  \item{maturity}{ maturity of the yield curve of which want to return the interest rates.}
}
\details{
\code{Coeff} is a vector or matrix of the four coefficients of the Nelson-Siegel's model: \eqn{(\beta_0; \beta_1; \beta_2; \lambda)}.
}
\value{
Return interest rates in matrix object with number of rows equal to \code{nrow(betaCoeff)} and number of columns equal to \code{length(maturity)}.
}
\references{ 
Diebold, F.X. and Li, C. (2006), Forecasting the Term Structure of Government Bond Yields, \emph{Journal of Econometrics}, \bold{130}, 337-364. 

Diebold, F.X., Ji, L. and Li, C. (2006), A Three-Factor Yield Curve Model: Non-Affine Structure, Systematic Risk Sources, and Generalized Duration, in L.R. Klein (ed.), \emph{Long-Run Growth and Short-Run Stabilization: Essays in Memory of Albert Ando}. Cheltenham, U.K.: Edward Elgar, 240-274.

Nelson, C.R., and A.F. Siegel (1987), Parsimonious Modeling of Yield Curve, \emph{The Journal of Business}, \bold{60}, 473-489.
}
\author{ Sergio Salvino Guirreri}

\examples{ 
data(FedYieldCurve)
maturity.Fed <- c(3/12, 0.5, 1,2,3,5,7,10)
NSParameters <- Nelson.Siegel( rate = first(FedYieldCurve,'10 month'), maturity=maturity.Fed )
y <- NSrates(NSParameters[5,],maturity.Fed)
plot(maturity.Fed,FedYieldCurve[10,],main="Fitting Nelson-Siegel yield curve", type="o")
lines(maturity.Fed,y, col=2)
legend("topleft",legend=c("observed yield curve","fitted yield curve"),
col=c(1,2),lty=1)
grid()
}
