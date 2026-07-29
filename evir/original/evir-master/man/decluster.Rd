\name{decluster}
\alias{decluster}
\title{Decluster Point Process}
\description{
  Declusters clustered point process data so that Poisson
  assumption is more tenable over a high threshold.
}
\usage{
decluster(series, run = NA, picture = TRUE)
}
\arguments{
  \item{series}{a numeric vector of threshold exceedances
    with a \code{times} attribute which should be a numeric
    vector containing either the indices or the times/dates
    of each exceedance (if times/dates, the attribute should
    be an object of class \code{"POSIXct"} or an object that
    can be converted to that class; see
    \code{\link{as.POSIXct}})}
  \item{run}{parameter to be used in the runs method; any two
    consecutive threshold exceedances separated by more than this
    number of observations/days are considered to belong to
    different clusters}  
  \item{picture}{whether or not a picture of declustering should
    be drawn}
}
\value{
  The declustered object. 
}
\section{References}{
  Embrechts, P., Klueppelberg, C., Mikosch, T. (1997).
  \emph{Modelling Extremal Events}. Springer.
  Chapter 8, 413--429. 
}
\seealso{\code{\link{pot}}, \code{\link{exindex}},
  \code{\link{as.POSIXct}}}

\examples{
# decluster the 200 exceedances of a particular threshold in 
# the negative BMW data
data(bmw)
out <- pot(-bmw, ne = 200) 
decluster(out$data, 30) 
}
\keyword{manip}
