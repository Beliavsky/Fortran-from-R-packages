\name{ECBYieldCurve}
\alias{ECBYieldCurve}
\docType{data}
\title{ Yield curve data spot rate, AAA-rated bonds, maturities from 3 months to 30 years}
\description{
Government bond, nominal, all triple A issuer companies. The maturities are 3 and 6 months and from 1 year to 30 years with frequency business day, provided by European Central Bank.
The range date is from 2006-12-29 to 2009-07-24.
}
\usage{data(ECBYieldCurve)}
\format{
It is an \code{xts} object with 32 interest rate at different maturities and 655 observations.
}

\source{
\acronym{ECB}: \url{https://www.ecb.europa.eu/stats/financial_markets_and_interest_rates/euro_area_yield_curves/html/index.en.html}.\cr
}

\examples{
### plot ECB Yield Curve ###
data(ECBYieldCurve)

first(ECBYieldCurve,'3 day')
last(ECBYieldCurve,'3 day')

mat.ECB <- tau <- c(3/12,6/12,1:30)

par(mfrow=c(2,3))
for( i in c(1,2,3,653,654,655) ){
plot(mat.ECB, ECBYieldCurve[i,], type="o", xlab="Maturity in years", ylab="IR values")
title(main=paste("European Central Bank yield curve obeserved at",time(ECBYieldCurve[i], sep=" ") ))
grid()
}
}
\keyword{datasets}
