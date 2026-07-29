\name{plotHistAcd}
\alias{plotHistAcd}
\title{
Mean duration plot
}
\description{
Plots the mean duration over time at chosen interval length
}
\usage{
plotHistAcd(durations, windowunit = "mins", window = 1)
}

\arguments{
  \item{durations}{
a data.frame containing the durations and their time stamps.
}
  \item{windowunit}{
the unit of the time interval. One of "secs", "mins", "hours", or "days".
}
  \item{window}{
a positive integer giving the length of the interval.
}
}
\author{
Markus Belfrage
}
\examples{
data(durData)
plotHistAcd(durData, windowunit = "days", window = 1)

\dontrun{

plotHistAcd(durData, windowunit = "mins", window = 30)
}
}
