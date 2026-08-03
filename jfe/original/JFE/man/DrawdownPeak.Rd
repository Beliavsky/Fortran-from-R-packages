\name{DrawdownPeak}
\alias{DrawdownPeak}

\encoding{latin1}

\title{Drawdawn peak of the return distribution}
\usage{
DrawdownPeak(R)
}
\arguments{
\item{R}{an xts, vector, matrix, data frame, timeSeries or zoo object of
asset returns}

}
\description{
Drawdawn peak is for each return its drawdown since the previous peak
}

\examples{
  data(assetReturns)
	R=assetReturns[, -29]
# Not run
# DrawdownPeak(R)


}

\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}
