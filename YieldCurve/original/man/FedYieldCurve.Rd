\name{FedYieldCurve}
\alias{FedYieldCurve}
\docType{data}
\title{Federal Reserve interest rates }
\description{
The data-set contains the interest rates of the Federal Reserve, from January 1982 to December 2012.
The interest rates are Market yield on U.S. Treasury securities constant maturity (CMT)
(more information on the Treasury yield curve can be found at the following website
\url{https://home.treasury.gov/policy-issues/financing-the-government/interest-rate-statistics})
at different maturities (3 months, 6 months, 1 year, 2 years, 3 years, 5 years, 7 years and 10 years),
quoted on investment basis and have been gathered with monthly frequency.
}
\usage{data(FedYieldCurve)}
\format{
  An object with class attributes \code{xts}.
}
\source{
\acronym{FED}: \url{https://www.federalreserve.gov/datadownload/Build.aspx?rel=H15}.
}

\examples{
require(xts)
require(YieldCurve)
data(FedYieldCurve)

first(FedYieldCurve,'3 month')
last(FedYieldCurve,'3 month')
mat<-c(3/12, 0.5, 1,2,3,5,7,10)

par(mfrow=c(2,3))
for( i in c(1,2,3,370,371,372) ){
plot(mat, FedYieldCurve[i,], type="o", xlab="Maturity in years", ylab="IR values")
title(main=paste("Federal Reserve yield curve obeserved at",time(FedYieldCurve[i], sep=" ") ))
grid()
}
}
\keyword{datasets}