\name{maxDrawdown}
\alias{maxDrawdown}
\encoding{latin1}

\title{caclulate the maximum drawdown from peak equity}
\usage{
maxDrawdown(R, geometric = TRUE, invert = TRUE)
}
\arguments{
\item{R}{an xts, vector, matrix, data frame, timeSeries or zoo object of
asset returns}

\item{geometric}{utilize geometric chaining (TRUE) or simple/arithmetic chaining (FALSE) to aggregate returns,
default TRUE}

\item{invert}{TRUE/FALSE whether to invert the drawdown measure.  see
Details.}


}
\description{
To find the maximum drawdown in a return series, we need to first calculate
the cumulative returns and the maximum cumulative return to that point.  Any
time the cumulative returns dips below the maximum cumulative returns, it's
a drawdown.  Drawdowns are measured as a percentage of that maximum
cumulative return, in effect, measured from peak equity.
}
\details{
The option to \code{invert} the measure should appease both academics and
practitioners. The default option \code{invert=TRUE} will provide the
drawdown as a positive number.  This should be useful for optimization
(which usually seeks to minimize a value), and for tables (where having
negative signs in front of every number may be considered clutter).
Practitioners will argue that drawdowns denote losses, and should be
internally consistent with the quantile (a negative number), for which
\code{invert=FALSE} will provide the value they expect.  Individually,
different preferences may apply for clarity and compactness.  As such, we
provide the option, but make no value judgment on which approach is
preferable.
}
\examples{

  data(assetReturns)
	R=assetReturns[, -29]

  maxDrawdown(R)

}
\references{
Bacon, C. \emph{Practical Portfolio Performance Measurement and
Attribution}. Wiley. 2004. p. 88 \cr
See also package \code{PerformanceAnalytics}.
}

\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}
