\name{plotDescTrans}
\alias{plotDescTrans}

\title{
Transactions plots 
}
\description{
Plots (1) the price over time, (2) volume traded over time for a given interval, and (3) number of transactions over time for a given interval. 
}
\usage{
plotDescTrans(trans, windowunit = "hours", window = 1)
}
\arguments{
  \item{trans}{
a data.frame with the column 'time', 'price', and 'volume'. Currently only works if all of those are available.
}
  \item{windowunit}{
the unit of the time interval. One of "secs", "mins", "hours", or "days".
}
  \item{window}{
a positive integer giving the length of the interval.
}
}
\examples{
\dontrun{
data(transData)
plotDescTrans(transData, windowunit = "hours", window = 1)}
}