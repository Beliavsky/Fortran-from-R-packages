\name{testSTACD}
\alias{testSTACD}

\title{
LM test against Smooth Transition ACD models (Meitz and Terasvirta, 2006)
}
\description{
Tests if the alpha parameters and the constant should be varying with the value of the lagged durations, according to a logistic transition function (currently only works for \code{model = "ACD"}).
}
\usage{
testSTACD(fitModel, K = 2, robust = TRUE)
}
\arguments{
  \item{fitModel}{
a fitted ACD model, i.e. an object of class "acdFit".
}
  \item{K}{
the order of the logistic transition function used for the alternative hypothesis.
}
  \item{robust}{
if \code{TRUE} the LM statistic will be calculated using the "robust" version, making its asymptotic behavior unaffected by possible misspecification of the error term distribution (Meitz and Terasvirta, 2006).
}
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

\seealso{
\code{\link{testRmACD}}, \code{\link{testTVACD}}.
}
\examples{
data(adjDurData)
fitModel3000obs <- acdFit(adjDurData[1:3000,])
testSTACD(fitModel3000obs, K = 2, robust = TRUE)

}

\references{
Meitz, M. and Terasvirta, T. (2006) Evaluating models of autoregressive conditional duration. Journal of Business and Economic Statistics 24: 104-124.
}