\name{exindex}
\alias{exindex}
\title{Estimate Extremal Index}
\description{
  Plot estimates of extremal index using the blocks method. 
}
\usage{
exindex(data, block, start = 5, end = NA, reverse = FALSE,
    auto.scale = TRUE, labels = TRUE, \dots)
}
\arguments{
  \item{data}{data vector (raw values not block maxima).}
  \item{block}{the block size. A numeric value is interpreted as
    the number of data values in each successive block.
    All the data is used, so the last block may not contain
    \code{block} observations.
    If the \code{data} has a \code{times} attribute containing (in
    an object of class \code{"POSIXct"}, or an object that can be
    converted to that class; see \code{\link{as.POSIXct}}) the
    times/dates of each observation, then \code{block} may instead
    take the character values \code{"month"}, \code{"quarter"},
    \code{"semester"} or \code{"year"}.}
  \item{start}{lowest value of K at which to plot a point; K is the
    number of blocks in which a specified threshold is exceeded}
  \item{end}{highest value of K at which to plot a point}
  \item{reverse}{whether plot is to be by increasing threshold
    (\code{TRUE}) or increasing K value (\code{FALSE})}
  \item{auto.scale}{whether or not plot should be automatically
    scaled; if not, \code{xlim} and \code{ylim} graphical parameters
    may be entered}
  \item{labels}{whether or not axes should be labelled}
  \item{\dots}{other graphics parameters}
}
\value{
  A table of results is returned invisibly.
}
\references{
  Embrechts, P., Klueppelberg, C., Mikosch, T. (1997).
  \emph{Modelling Extremal Events}. Springer. Chapter 8, 413-429.  
}
\seealso{\code{\link{gev}}, \code{\link{hill}},
  \code{\link{as.POSIXct}}}

\examples{
\dontrun{data(bmw)}
\dontrun{exindex(bmw, 100)} 
\dontrun{exindex(-bmw, 100)} 
# calculate extremal index for the right and left tails of the BMW 
# log returns 
}
\keyword{hplot}
