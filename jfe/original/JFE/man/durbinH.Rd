\name{durbinH}
\alias{durbinH}
\encoding{latin1}
\title{calculate Sortino Ratio of performance over downside risk}
\usage{
durbinH(model,Ly.label="ar1")
}
\arguments{
\item{model}{Estimation result object by, for example, dynlm or arima}

\item{Ly.label}{Label of AR(1) dependent variable}

}
\description{
durbinH test the serial correlation of a dynamic regression, namely, a time
series regression with AR(1) terms. It modifies Durbin-Watson statictic.
}

\examples{

  data(macrodata)
  y=ts(macrodata[,"OECD"],start=c(1961,1),freq=12) #OECD business cycle index

 out1=dynlm::dynlm(y~L(y,1)+season(y)+trend(y))
 durbinH(out1,Ly.label="L(y, 1)")

}

\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

