\name{acf_acd}
\alias{acf_acd}
\title{
Autocorrelation function plots for ACD models
}
\description{
plots the ACF (Auto Correlation Function) for the durations, diurnally adjusted durations, and residuals.
}
\usage{
acf_acd(fitModel = NULL, conf_level = 0.95, max = 50, min = 1,
        plotDur = TRUE, plotAdjDur = TRUE, plotResi = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitModel}{
a fitted model of class "acdFit", or a data.frame containing at least one the columns "durations", "adjDur", or "residuals". Can also be a vector of durations or residuals.
}
  \item{conf_level}{
the confidence level of the confidence bands
}
  \item{max}{
the largest lag to plot
}
  \item{min}{
the smallest lag to plot
}
  \item{plotDur, plotAdjDur, plotResi}{
logical falgs. If \code{FALSE}, the respective ACF wont be plotted.
}
}
\value{
returns a data.frame with the values of the sample autocorrelations for each lag and variable.
}
\author{
Markus Belfrage
}
\examples{

data(adjDurData)
fitModel <- acdFit(adjDurData)
acf_acd(fitModel, conf_level = 0.95, max = 50, min = 1)

data(durData)
f <- acf_acd(durData)
f

}