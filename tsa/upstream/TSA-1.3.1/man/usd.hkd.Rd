\name{usd.hkd}
\alias{usd.hkd}
\docType{data}
\title{Daily US Dollar to Hong Kong Dollar Exchange Rates}
\description{
Daily USD/HKD
(US dollar to Hong Kong dollar) exchange rate from January 1, 2005 to March 7, 2006
}
\usage{data(usd.hkd)}
\format{
  A data frame with 431 observations on the following 6 variables.
  \describe{
    \item{\code{r}}{daily returns of USD/HKD exchange rates}
    \item{\code{v}}{estimated conditional variances based on an AR(1)+GARCH(3,1) model}
    \item{\code{hkrate}}{daily USD/HKD exchange rates}
    \item{\code{outlier1}}{dummy variable of day 203, corresponding to July 22, 2005}
    \item{\code{outlier2}}{dummy variable of day 290, another possible outlier}
    \item{\code{day}}{calendar day}
  }
}
\source{
\code{http://www.oanda.com/convert/fxhistory}
}
\references{ "Time Series Analysis, with Applications in R" by J.D. Cryer and K.S. Chan}
\examples{
data(usd.hkd)
## maybe str(usd.hkd) ; plot(usd.hkd) ...
}
\keyword{datasets}
