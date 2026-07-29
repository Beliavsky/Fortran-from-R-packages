\name{qqplotAcd}
\alias{qqplotAcd}

\title{

Quantile-Quantile plot of the residuals
}
\description{
Plots a QQ-plot of the residuals and the theoretical quantiles implied by the model estimates.
}
\usage{
qqplotAcd(fitModel, xlim = NULL, ylim = NULL)
}
\arguments{
  \item{fitModel}{
a fitted ACD model, i.e. an object of class "acdFit"
}
  \item{xlim}{
an optional vector of limits for the x-axis
}
  \item{ylim}{
an optional vector of limits for the y-axis
}
}
\examples{
data(adjDurData)
fitModelExp <- acdFit(adjDurData, dist = "exp")
qqplotAcd(fitModelExp)
}