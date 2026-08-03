\name{AppraisalRatio}
\alias{AppraisalRatio}
\title{Appraisal ratio of the return distribution}

\encoding{latin1}

\usage{
AppraisalRatio(Ra, Rb, Rf = 0, method = c("appraisal", "modified",
  "alternative"))
}
\arguments{
\item{Ra}{an xts, vector, matrix, data frame, timeSeries or zoo object of
asset returns}

\item{Rb}{return vector of the benchmark asset}

\item{Rf}{risk free rate, in same period as your returns}

\item{method}{is one of "appraisal" to calculate appraisal ratio, "modified" to
calculate modified Jensen's alpha or "alternative" to calculate alternative
Jensen's alpha.}

}
\description{
Appraisal ratio is the Jensen's alpha adjusted for specific risk. The numerator
is divided by specific risk instead of total risk.
}
\details{
Modified Jensen's alpha is Jensen's alpha divided by beta.

Alternative Jensen's alpha is Jensen's alpha divided by systematic risk.

\deqn{Appraisal ratio = \frac{\alpha}{\sigma_{\epsilon}}}{Appraisal ratio = Jensen's alpha / specific risk}

\deqn{Modified Jensen's alpha = \frac{\alpha}{\beta}}{Modified Jensen's alpha = Jensen's alpha / beta}

\deqn{Alternative Jensen's alpha = \frac{\alpha}{\sigma_S}}{Alternative Jensen's alpha = Jensen's alpha / systematic risk}

where \eqn{alpha} is the Jensen's alpha, \eqn{\sigma_{epsilon}} is the specific risk,
\eqn{\sigma_S} is the systematic risk.
}
\examples{

    data(assetReturns)
	Ra=assetReturns[, -29]
	Rb=assetReturns[,29] #DJI
AppraisalRatio(Ra, Rb, method="appraisal")

}
\references{
Carl Bacon, \emph{Practical portfolio performance measurement
and attribution}, second edition 2008 p.77.\cr
See also package \code{PerformanceAnalytics}.
}
\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}
