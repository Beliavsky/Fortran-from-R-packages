\name{resiDensityAcd}
\alias{resiDensityAcd}

\title{
Residual Density Histogram
}
\description{
Plots a density histogram of the residuals and superimposes the density implied by the model estimates.
}
\usage{
resiDensityAcd(fitModel, xlim = NULL, binwidth = .1, density = FALSE)
}
\arguments{
  \item{fitModel}{
a fitted ACD model, i.e. an object of class "acdFit"
}
  \item{xlim}{
an optional vector of limits for the x-axis
}
  \item{binwidth}{
the width of the bins of the density histogram.
}
  \item{density}{
if \code{TRUE} a kernel density estimate will be added
}
}
\author{
Markus Belfrage
}
\examples{
\dontrun{
data(adjDurData)
fitModelBurr <- acdFit(adjDurData, dist = "burr")
resiDensityAcd(fitModelBurr)}
}