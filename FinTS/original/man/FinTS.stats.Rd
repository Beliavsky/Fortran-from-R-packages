\name{FinTS.stats}
\alias{FinTS.stats}
\title{ Financial Time Series summary statistics }
\description{
  Summary statistics as in Table 1.2, Tsay (2005), including the start
  date, number of observations, mean, standard deviation, skewness,
  excess kurtosis, min and max.  
}
\usage{
FinTS.stats(x) 
}
\arguments{
  \item{x}{
    A univariate object of class 'zoo' 
  }
}
\references{
  Ruey Tsay (2005)
  Analysis of Financial Time Series, 2nd ed. (Wiley, p. 11)
}
\seealso{
%  \code{\link{Table1.2}}
  \code{\link[zoo]{index}},
  \code{\link{sum}},
  \code{\link{is.na}},
  \code{\link{mean}},
  \code{\link{sd}},
  \code{\link[e1071]{skewness}},
  \code{\link[e1071]{kurtosis}},
  \code{\link{min}},
  \code{\link{max}}
}
\examples{
FinTS.stats(rep(1, 5))

data(d.c8603)
FinTS.stats(100*d.c8603[, "C"])

## The following generates an error, because FinTS.stats
## expects a vector of class 'zoo', and d.c8603 is a matrix
#FinTS.stats(100*d.c8603) 
}
\keyword{univar}
