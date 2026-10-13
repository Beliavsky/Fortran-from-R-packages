\name{getParamsBaseScenarios}
\alias{getParamsBaseScenarios}
\alias{getParamsBaseScenarios,Scenarios-method}
\title{getParamsBaseScenarios method}
\description{
  Get a list containing the horizon and number of scenarios
  for a Scenarios object.
}
\examples{
scenarios1 <- new("Scenarios")
scenarios1 <- setParamsBaseScenarios(scenarios1, horizon=5, nScenarios=10)
getParamsBaseScenarios(scenarios1)
}

