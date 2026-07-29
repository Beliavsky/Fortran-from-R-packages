\name{ciTarThd}
\alias{ciTarThd}
\alias{print.ciTarThd}
\alias{plot.ciTarThd}
\title{
Threshold Selection for Threshold Cointegration Regression
}
\description{
Select the best threshold for threshold cointegration regression by sum of squared errors
}
\usage{
  ciTarThd(y, x, model = c('tar', 'mtar'), lag = 1, 
  th.range = 0.15, digits = 3)
}
\arguments{
  \item{y}{dependent or left-side variable for the long-run regression.}
  \item{x}{independent or right-side variable for the long-run regression.}
  \item{model}{a choice of two models, either tar or mtar.}
  \item{lag}{number of lags.} 
  \item{th.range}{the percentage of observations to be excluded from the search.}
  \item{digits}{number of digits used in rounding outputs.}
}
\details{
The best threshold is determined by fitting the regression for possible threshold values, sorting the results by sum of squared errors (SSE), and selecting the best with the lowest SSE. To have sufficient observations on either side of the threshold value, certain percentage of observations on the top and bottoms are excluded from the search path. This is usually set as 0.15 by the \code{th.range} (Chan 1993).
}
\value{
Return a list object of class "ciTarThd" with the following components: 
  \item{model}{model choice} 
  \item{lag}{number of lags} 
  \item{th.range}{the percentage of observations excluded} 
  \item{th.final}{the best threshold value}
  \item{ssef}{the best (i.e., lowest) value of SSE}
  \item{obs.tot}{total number of observations in the raw data}
  \item{obs.CI}{number of observations used in the threshold cointegration regression}      
  \item{basic}{a brief summary of the major outputs}
  \item{path}{a data frame of the search record (number of regression, threshold value, SSE, AIC, and BIC values).}
}

\section{Methods}{
  Two methods are defined as follows:
  \describe{
      \item{\code{print}:}{ This shows the \code{basic} component in the returned list object.}

    \item{\code{plot}:}{ plotting three graphs in one window; they reveals the relationship between SSE (sum of squared errors), AIC, BIC and the threshold values. The best threshold value is associated with the lowest SSE value.}
  }
} 

\references{
Chan, K.S. 1993. Consistency and limiting distribution of the least squares estimator of a threshold autoregressive model. The Annals of Statistics 21(1):520-533.

Enders, W., and C.W.J. Granger. 1998. Unit-root tests and asymmetric adjustment with an example using the term structure of interest rates. Journal of Business & Economic Statistics 16(3):304-311.
}

\author{Changyou Sun (\email{edwinsun258@gmail.com})}
\seealso{
\code{\link{ciTarFit}}; and \code{\link{ciTarLag}} 
}
\examples{
# see example at daVich
}
\keyword{regression}