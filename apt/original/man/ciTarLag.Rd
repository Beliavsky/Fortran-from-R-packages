\name{ciTarLag}
\alias{ciTarLag}
\alias{print.ciTarLag}
\alias{plot.ciTarLag}
\title{
  Lag Selection for Threshold Cointegration Regression
}
\description{
  Select the best lag for threshold cointegration regression by AIC and BIC
}
\usage{
  ciTarLag(y, x, model = c("tar", "mtar"), maxlag = 4, 
  thresh = 0, adjust = TRUE) 
}
\arguments{
  \item{y}{dependent or left-side variable for the long-run regression.}
  \item{x}{independent or right-side variable for the long-run regression.}
  \item{model}{a choice of two models, either tar or mtar.}
  \item{maxlag}{maximum number of lags allowed in the search process.}
  \item{thresh}{a threshold value.}
  \item{adjust}{logical value (default of TRUE) of whether to adjust the window widths so all regressions by lag have the same number of observations}
}

\details{
Estimate the threshold cointegration regressions by lag and then select the best regression by AIC or BIC value. The longer the lag, the smaller the number of observations availabe for estimation. If the windows of regressions by lag are not ajusted, the maximum lag is usually the best lag by AIC or BIC. Theorectially, AIC and BIC from different models should be compared on the basis of the same observation numbers (Ender 2004). \code{adjust} shows the effect of this adjustment on the estimation window. By default, the value of \code{adjust} should be TRUE. 
}
\value{
Return a list object of class "ciTarLag" with the following components: 
  \item{path}{a data frame of model criterion values by lag, including \code{lag} for the current lag, \code{totObs} for total observations in the raw data, \code{coinObs} for observations used in the cointegration regression, \code{sse} for the sum of squared errors, \code{aic} for AIC value, \code{bic} for BIC value, \code{LB4} for the p-value of Ljung_Box Q statistic with 4 autocorrelation coefficients, \code{LB8} with 8 coefficients, \code{LB12} for Q statistic with 12 coefficients}
   
  \item{out}{a data frame of the final model selection, including the values of model, maximum lag, threshold value, best lag by AIC, best lag by BIC} 
}

\section{Methods}{
  Two methods are defined as follows:
  \describe{
    \item{\code{print}:}{ This shows the \code{out} component in the returned list.}
    \item{\code{plot}:}{ This demonstrates the trend of AIC and BIC changes of threshold cointegration regressions by lag. It facilitates the selection of the best lag for a threshold cointegration model.}
  }
} 

\references{
Enders, W. 2004. Applied Econometric Time Series. John Wiley & Sons, Inc., New York. 480 P.

Enders, W., and C.W.J. Granger. 1998. Unit-root tests and asymmetric adjustment with an example using the term structure of interest rates. Journal of Business & Economic Statistics 16(3):304-311.
}
\author{Changyou Sun (\email{edwinsun258@gmail.com})}
\seealso{
\code{\link{ciTarFit}}; and \code{\link{ciTarThd}}; 
}
\examples{
# see example at daVich
}
\keyword{regression}