\name{diurnalAdj}
\alias{diurnalAdj}
\title{
Dirunal adjustment for durations
}
\description{
Performs a diurnal adjustment of the durations, i.e. removes a daily seasonal component. Four different methods of diurnal adjustment are available, namely "cubicSpline", "supsmu" (Friedman's SuperSmoother), "smoothSpline" (smoothed version of the cubic spline), or "FFF" (Flexible Fourier Form).
}
\usage{
diurnalAdj(dur, method = "cubicSpline", nodes = c(seq(600, 1105, 60), 1105),
aggregation = "all", span = "cv", spar = 0, Q = 4, returnSplineFnc = FALSE)
}
\arguments{
  \item{dur}{
a data.frame containing the columns \code{durations}, containing durations, and \code{time}, containing the time stamps.
}
  \item{method}{
the method used. One of "cubicSpline", "supsmu", "smoothSpline", or "FFF".
}
  \item{nodes}{
only for \code{method = "cubicSpline"} or \code{method = "smoothSpline"}. A vector of nodes to use for the spline function, in the unit minutes after midnight. The first and last element of the vector must be the start and end of the trading day. The nodes given are actually the limits of intervalls, of wich the midpoints will be set as the nodes using the means of the intervals. 
}
  \item{aggregation}{
what type of aggregation to use. Either "weekdays", "all", or "none". If for example "weekdays" is chosen, all Mondays will have the same daily seasonal component, and so on. 
}
  \item{span}{
argument passed to supsmu if \code{method = "supsmu"} were chosen. Affects the smoothness of the curve, see \link{supsmu}.
}
  \item{spar}{
argument passed to smooth.spline if \code{method = "smooth.spline"} were chosen. Affects the smoothness of the curve, see \link{smooth.spline}.
}
  \item{Q}{
number of trigonometric function pairs for \code{method = "FFF"}. 
}
  \item{returnSplineFnc}{
if \code{TRUE} instead or returning the adjusted durations a list of spline objects will be returned, containing the coefficents of the spline function. Only available for \code{method = "cubicSpline"}. 
}
}
\value{
\item{if \code{returnSplineFnc = FALSE} (default)}{
the input data.frame 'dur' with an added column of the diurnally adjusted durations called 'adjDur'.
}
\item{if \code{returnSplineFnc = TRUE}}{
a list of spline objects containing the coefficents of the spline function.
}
}
\author{
Markus Belfrage
}
\examples{
data(durData)
diurnalAdj(durData, aggregation = "none", method = "supsmu")

\dontrun{

head(durData)
f <- diurnalAdj(durData, aggregation = "weekdays", method = "FFF", Q = 3)
head(f)

f <- diurnalAdj(durData, aggregation = "all", returnSplineFnc = TRUE)
f
}
}
