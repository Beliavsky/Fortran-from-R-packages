\name{setParamsBaseScenarios}
\alias{setParamsBaseScenarios}
\alias{setParamsBaseScenarios,Scenarios-method}
\title{setParamsBaseScenarios method}
\arguments{
  \item{horizon}{Horizon for the projection (in years)}

  \item{nScenarios}{Number of scenarios}
}
\description{
  Set the horizon and nScenarios parameters of the
  [\code{\linkS4class{ParamsScenarios}}] sub-object of a
  Scenarios object
}
\examples{
scenarios1 <- new("Scenarios")
scenarios1 <- setParamsBaseScenarios(scenarios1, horizon=5, nScenarios=10)
}

