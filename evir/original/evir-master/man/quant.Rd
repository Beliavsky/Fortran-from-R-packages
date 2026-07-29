\name{quant}
\alias{quant}
\title{Plot of GPD Tail Estimate of a High Quantile}
\description{
  Creates a plot showing how the estimate of a high quantile in
  the tail of a dataset based on the GPD approximation varies with
  threshold or number of extremes. 
}
\usage{
quant(data, p = 0.99, models = 30, start = 15, end = 500, reverse =
    TRUE, ci = 0.95, auto.scale = TRUE, labels = TRUE, \dots)
}
\arguments{
  \item{data}{numeric vector of data}
  \item{p}{desired probability for quantile estimate (e.g. 0.99
    gives 99th percentile)}
  \item{models}{number of consecutive gpd models to be fitted}
  \item{start}{lowest number of exceedances to be considered}
  \item{end}{maximum number of exceedances to be considered}
  \item{reverse}{should plot be by increasing threshold
    (\code{TRUE}) or number of extremes (\code{FALSE})}
  \item{ci}{probability for asymptotic confidence band; for no
    confidence band set to zero}
  \item{auto.scale}{whether or not plot should be automatically scaled;
    if not, xlim and ylim graphical parameters may be entered}
  \item{labels}{whether or not axes should be labelled}
  \item{\dots}{other graphics parameters}
}
\value{
  A table of results is returned invisibly.
}
\details{
  For every model \code{gpd} is called. Evaluation may be slow. 
  Confidence intervals by the Wald method (which is fastest). 
}
\seealso{\code{\link{gpd}}, \code{\link{plot.gpd}},
  \code{\link{gpd.q}}, \code{\link{shape}}}

\examples{
\dontrun{data(danish)}
\dontrun{quant(danish, 0.999)} 
# Estimates of the 99.9th percentile of the Danish losses using 
# the GPD model with various thresholds 
}
\keyword{hplot}
