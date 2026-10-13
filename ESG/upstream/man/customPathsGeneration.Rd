\name{customPathsGeneration}
\alias{customPathsGeneration}
\alias{customPathsGeneration,Scenarios-method}
\title{customPathsGeneration method}
\arguments{
  \item{type}{The name of the asset for which a projection
  has to be proceeded. Can be 'shortRate', 'stock',
  'realEstate', 'liquiditySpread' or 'defaultSpread'. If
  NULL, all assets will be projected.}
}
\description{
  Proceed to the projection using the parameters that were
  previously set into the Scenarios objet.
}
\examples{
objScenario  <- new("Scenarios")
# Basic scenario's parameters setting
objScenario  <- setParamsBaseScenarios(objScenario, horizon = 10, nScenarios = 1000)
# Risk factors parameters setting
objScenario  <- setRiskParamsScenariosrt(objScenario, vol = .1, k = 2)
objScenario  <- setRiskParamsScenariosS(objScenario, vol = .1, k = 2,
volStock = .2, stock0 = 100, rho=.5)
objScenario  <- setRiskParamsScenariosliqSpr(objScenario, eta=.05, liquiditySpread0=.01)
objScenario  <- setRiskParamsScenariosdefSpr(objScenario, volDefault=.2,
defaultSpread0=.01, alpha=.1, beta=1)
# Forward and ZC rates setting
data(ZC)
objScenario  <- setForwardRates(objScenario, ZC, horizon=10)
objScenario  <- setZCRates(objScenario, ZC, horizon=10)
# Projection
objScenario  <- customPathsGeneration(objScenario, type="shortRate")
objScenario  <- customPathsGeneration(objScenario, type="stock")
objScenario  <- customPathsGeneration(objScenario, type="defaultSpread")
objScenario  <- customPathsGeneration(objScenario, type="liquiditySpread")
}

