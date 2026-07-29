\name{testTVACD}
\alias{testTVACD}

\title{
LM test against Time-Varying ACD models (Meitz and Terasvirta, 2006)
}
\description{
Tests if the parameters are time-varying (currently only works for \code{model = "ACD"}).
}
\usage{
testTVACD(fitModel, K = 2, type = "total", robust = TRUE)
}
\arguments{
  \item{fitModel}{
a fitted ACD model, i.e. an object of class "acdFit".
}
  \item{K}{
the order of the logistic transition function used for the alternative hypothesis.
}
  \item{type}{
either "total" or "intraday". If "total", the possible time varying parameters under the alternative varies over the total time of the sample, whereas for "intraday", the time variable is time of the day. See 'Details'
}
  \item{robust}{
if \code{TRUE} the LM statistic will be calculated using the "robust" version, making its asymptotic behavior unaffected by possible misspecification of the error term distribution (Meitz and Terasvirta, 2006).
}
}
\details{
This function tests the fitted standard ACD model against the TVACD model of Meitz and Terasvirta (2006). The TVACD model lets the ACD parameters vary over time by a logistic transition function.

In one specification, the time variable is total time, and a test rejecting the null in favor of this alternative specification would indicate that the ACD parameters are changing over time over the total sample.

The other specification lets the parameters be intraday varying, by letting the transition variable be the time of the day. Failing this test could indicate that the diurnal adjustment was inadequate at removing any diurnal component. 
}
\value{
a list of:
\item{chi2}{
the value of the LM statistic.
}
\item{pv}{
the pvalue of the test statistic.
}
}
\references{
Meitz, M. and Terasvirta, T. (2006).
\emph{Evaluating models of autoregressive conditional duration}. Journal of Business and Economic Statistics 24: 104-124. 
}
\author{
Markus Belfrage
}
\seealso{
\code{\link{testRmACD}}, \code{\link{testSTACD}}.
}
\examples{
data(adjDurData)
fitModel5000obs <- acdFit(adjDurData[1:5000,])
testTVACD(fitModel5000obs, K = 2, type = "total", robust = TRUE)

testTVACD(fitModel5000obs, K = 2, type = "intraday", robust = TRUE)

}