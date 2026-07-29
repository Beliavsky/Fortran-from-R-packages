\name{shape}
\alias{shape}
\title{Plot for GPD Shape Parameter}
\description{
  Creates a plot showing how the estimate of shape varies
  with threshold or number of extremes. 
}
\usage{
shape(data, models = 30, start = 15, end = 500, reverse = TRUE,
    ci = 0.95, auto.scale = TRUE, labels = TRUE, \dots)
}
\arguments{
  \item{data}{numeric vector of data}
  \item{models}{number of consecutive gpd models to be fitted}
  \item{start}{lowest number of exceedances to be considered}
  \item{end}{maximum number of exceedances to be considered}
  \item{reverse}{should plot be by increasing threshold
    (\code{TRUE}) or number of extremes (\code{FALSE})}
  \item{ci}{probability for asymptotic confidence band; for no
    confidence band set to zero}
  \item{auto.scale}{whether or not plot should be automatically
    scaled; if not, \code{xlim} and \code{ylim} graphical
    parameters may be entered}
  \item{labels}{whether or not axes should be labelled}
  \item{\dots}{other graphics parameters}
}
\value{
  A table of results is returned invisibly.
}
\details{
  For every model \code{gpd} is called. Evaluation may be slow. 
}
\seealso{\code{\link{gpd}}, \code{\link{plot.gpd}},
  \code{\link{hill}}}

\examples{
\dontrun{data(danish)}
\dontrun{shape(danish)}
# Shape plot of heavy-tailed Danish fire insurance data 
}
\keyword{hplot}
