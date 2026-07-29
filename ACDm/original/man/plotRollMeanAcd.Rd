\name{plotRollMeanAcd}
\alias{plotRollMeanAcd}
\title{
Plots rolling means of durations
}
\description{
Plots rolling means of durations
}
\usage{
plotRollMeanAcd(durations, window = 500)
}
\arguments{
  \item{durations}{
a data.frame containing the column 'time' and 'durations'.
}
  \item{window}{
the length of the rolling window.
}
}
\examples{
data(durData)
plotRollMeanAcd(durData, window = 500)

}