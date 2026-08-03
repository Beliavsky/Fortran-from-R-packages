    \name{gas}
     \docType{data}
     \alias{gas}
     \title{Standardized residuals of weekly log-returns of gas and oil prices in Canada from 2008 to end of February 2011}
     \description{Data frame containg xres (standardized residuals of gas prices from a ARMA(2,2) model) and  yres (standardized residuals of oil prices from a ARMA(1,1)-GARCH(1,1) model).}
     \usage{data(gas)}
     \format{Residuals}
     \keyword{datasets}

\examples{
data(gas)
plot(gas$xres)
}
\keyword{datasets}

