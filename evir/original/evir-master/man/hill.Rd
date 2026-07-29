\name{hill}
\alias{hill}
\title{Create Hill Plot}
\description{
  Plot the Hill estimate of the tail index of heavy-tailed
  data, or of an associated quantile estimate.
}
\usage{
hill(data, option = c("alpha","xi","quantile"), start = 15,
    end = NA, reverse = FALSE, p = NA, ci = 0.95,
    auto.scale = TRUE, labels = TRUE, \dots)
}
\arguments{
  \item{data}{data vector}
  \item{option}{whether \code{"alpha"}, \code{"xi"} (1/alpha) or
    \code{"quantile"} (a quantile estimate) should be plotted}
  \item{start}{lowest number of order statistics at which to plot
    a point}
  \item{end}{highest number of order statistics at which to plot
    a point}
  \item{reverse}{whether plot is to be by increasing threshold
    (\code{TRUE}) or increasing number of order statistics
    (\code{FALSE})}
  \item{p}{probability required when option \code{"quantile"} is
    chosen}
  \item{ci}{probability for asymptotic confidence band; for no
    confidence band set \code{ci} to zero}
  \item{auto.scale}{whether or not plot should be automatically
    scaled; if not, \code{xlim} and \code{ylim} graphical
    parameters may be entered}
  \item{labels}{whether or not axes should be labelled}
  \item{\dots}{other graphics parameters}
}
\details{
  This plot is usually calculated from the alpha perspective.
  For a generalized Pareto analysis of heavy-tailed data using
  the \code{gpd} function, it helps to plot the Hill estimates
  for xi. 
}
\seealso{\code{\link{shape}}, \code{\link{quant}}}

\examples{
\dontrun{data(danish)}
\dontrun{hill(danish)} 
# Hill plot of heavy-tailed Danish fire insurance data 
\dontrun{hill(danish, option = "quantile", end = 500, p = 0.999)} 
# Hill plot of estimated 0.999 quantile of Danish fire insurance data 
}
\keyword{hplot}
