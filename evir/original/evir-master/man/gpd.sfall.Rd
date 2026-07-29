\name{gpd.sfall}
\alias{gpd.sfall}
\title{Add Expected Shortfall Estimates to a GPD Plot}
\description{
  Calculates expected shortfall (tail conditional expectation) 
  estimates and confidence intervals for high quantiles above the
  threshold in a GPD analysis, and adds a graphical representation
  to an existing plot.  
}
\usage{
gpd.sfall(x, pp, ci.p = 0.95, like.num = 50)
}
\arguments{
  \item{x}{a list object returned by \code{plot.gpd} or
    \code{tailplot}}
  \item{pp}{the desired probability for expected shortfall estimate
    (e.g. 0.99 for the 99th percentile)}
  \item{ci.p}{probability for confidence interval (must be less
    than 0.999)}
  \item{like.num}{number of times to evaluate profile likelihood}
}
\details{
  Expected shortfall is the expected size of the loss, given that a
  particular quantile of the loss distribution is exceeded. 
  The GPD approximation in the tail is used to estimate expected
  shortfall. The likelihood is reparametrised  in terms of the unknown
  expected shortfall and profile likelihood arguments are used to
  construct a confidence interval. 
}
\seealso{\code{\link{gpd}}, \code{\link{plot.gpd}},
  \code{\link{tailplot}}, \code{\link{gpd.q}}}

\examples{
\dontrun{data(danish)}
\dontrun{out <- gpd(danish, 10)}
\dontrun{tp <- tailplot(out)}
\dontrun{gpd.q(tp, 0.999)} 
# Estimates 99.9th percentile of Danish fire losses 
\dontrun{gpd.sfall(tp, 0.999)} 
# Estimates associated expected shortfall for Danish fire losses 
}
\keyword{iplot}
