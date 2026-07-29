\name{gpd.q}
\alias{gpd.q}
\title{Add Quantile Estimates to plot.gpd}
\description{
  Calculates quantile estimates and confidence intervals for
  high quantiles above the threshold in a GPD analysis, and
  adds a graphical representation to an existing plot.
}
\usage{
gpd.q(x, pp, ci.type = c("likelihood", "wald"), ci.p = 0.95,
    like.num = 50)
}
\arguments{
  \item{x}{a list object returned by \code{plot.gpd} or
    \code{tailplot}}
  \item{pp}{the desired probability for quantile estimate (e.g. 0.99
    for the 99th percentile)}
  \item{ci.type}{method for calculating a confidence interval:
    \code{"likelihood"} or \code{"wald"}}
  \item{ci.p}{probability for confidence interval (must be less
    than 0.999)}
  \item{like.num}{number of times to evaluate profile likelihood}
}
\details{
  The GPD approximation in the tail is used to estimate quantile. 
  The \code{"wald"} method uses the observed Fisher information
  matrix to calculate confidence interval. The \code{"likelihood"}
  method reparametrizes the likelihood in terms of the unknown
  quantile and uses profile likelihood arguments to construct a
  confidence interval. 
}
\seealso{\code{\link{gpd}}, \code{\link{plot.gpd}},
  \code{\link{gpd.sfall}}, \code{\link{tailplot}}}
\examples{
\dontrun{data(danish)}
\dontrun{out <- gpd(danish, 10)} 
\dontrun{tp <- tailplot(out)} 
\dontrun{gpd.q(tp, 0.999)} 
# Estimates 99.9th percentile of Danish fire losses 
}
\keyword{iplot}
