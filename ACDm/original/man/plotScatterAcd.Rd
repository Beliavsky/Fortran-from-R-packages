\name{plotScatterAcd}
\alias{plotScatterAcd}
\title{
Scatter plot for ACD models
}
\description{
Function to help scatter plot different variables of a fitted ACD model and superimposes a smoothed conditional mean using ggplot2. Can be used to investigate the possible need for non-linear models and issues with the diurnal adjustment.
}
\usage{
plotScatterAcd(fitModel, x = "muHats", y = "residuals", xlag = 0, ylag = 0,
                           colour = NULL, xlim = NULL, ylim = NULL, alpha = 1/10,
                           smoothMethod = "auto")
}
\arguments{
  \item{fitModel}{
a fitted model of class "acdFit"
}
  \item{x}{
the variable used on the x-axis. One of "muHats", "residuals", "durations", "adjDur", "dayTime", "time", or "index".
}
  \item{y}{
the variable used on the y-axis. One of "muHats", "residuals", "durations", "adjDur", "dayTime", "time", or "index".
}
  \item{xlag}{
number of lags used for the variable shown on the x-axis.
}
  \item{ylag}{
number of lags used for the variable shown on the y-axis.
}
  \item{colour}{
a possible third variable to be represented with a colour scale. One of "muHats", "residuals", "durations", "adjDur", "dayTime", or "time".
}
  \item{xlim}{
a vector of the limits of the x-axis to possibly zoom in on a certain region.
}
  \item{ylim}{
a vector of the limits of the y-axis to possibly zoom in on a certain region.
}
  \item{alpha}{
\code{alpha} parameter passed to ggplot2. For large data sets many data points will overlap. The alpha parameter can make the points transparent, making it easier to distinguish the density of different region. Takes the value between 1 (opaque) and 0 (completely transparent).
}
  \item{smoothMethod}{
value passed as \code{smooth} argument to ggplot2. See \link[ggplot2]{stat_smooth}.
}
}
\author{
Markus Belfrage
}
\examples{
\dontrun{

data(adjDurData)
# The mean residuals are too small for small values of the estimated conditional 
# mean, suggesting a need for a different conditional mean model specification:
fitModel <- acdFit(adjDurData)
plotScatterAcd(fitModel, x = "muHats", y = "residuals")
}
}