\name{tailplot}
\alias{tailplot}
\title{Plot Tail Estimate From GPD Model}
\description{
  Interacts with the output of gpd to produce a plot of the tail
  of the underlying distribution of the data.
  This is one of the options of \code{plot.gpd}, but tailplot enables
  the user to bypass the menu of the former. 
}
\usage{
tailplot(x, optlog = NA, extend = 1.5, labels = TRUE, \dots)
}
\arguments{
  \item{x}{a \code{gpd} object}
  \item{optlog}{optional argument giving a particular choice
    of logarithmic axes: \code{"x"} x-axis only; \code{"y"}
    y-axis only; \code{"xy"} both axes; \code{""} neither axis.}
  \item{extend}{optional argument expressing how far x-axis should
    extend as a multiple of the largest data value. This argument
    must take values greater than 1 and is useful for showing
    estimated quantiles beyond data.}
  \item{labels}{optional argument specifying whether or not axes
    should be labelled}
  \item{\dots}{other graphics parameters}
}
\value{
  A list object containing details of the plot is returned invisibly.
  This object should be used as the first argument of \code{gpd.q}
  or \code{gpd.sfall} to add quantile estimates or expected shortfall
  estimates to the plot.
}

\seealso{\code{\link{gpd}}, \code{\link{plot.gpd}},
  \code{\link{shape}}, \code{\link{quant}}}

\examples{
data(danish)
out <- gpd(danish, 10) 
\dontrun{tailplot(out)} 
}
\keyword{hplot}
