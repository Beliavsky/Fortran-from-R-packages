\name{setRiskParamsScenarios}
\alias{setRiskParamsScenarios}
\alias{setRiskParamsScenarios,Scenarios-method}
\title{setRiskParamsScenarios method}
\arguments{
  \item{vol}{Volatility for rates in vasicek model}

  \item{k}{k for rates in vasicek model}

  \item{volStock}{Volatility for UL in Black & Scholes
  model}

  \item{volRealEstate}{Volatility for real estate in Black
  & Scholes model}

  \item{volDefault}{Volatility for LMN model}

  \item{alpha}{alpha for LMN model}

  \item{beta}{beta Volatility for LMN model}

  \item{eta}{eta Volatility for LMN model}

  \item{rho}{Correlation between stock and short rates}

  \item{stock0}{UL initial value}

  \item{realEstate0}{Real estate initial value}

  \item{liquiditySpread0}{Initial liquidity for LMN model}

  \item{defaultSpread0}{Initial default spread for LMN
  model}
}
\description{
  Set all the risk parameters of a Scenarios object
  (contained in a [\code{\linkS4class{ParamsScenarios}}]
  sub-object)
}
\examples{
scenarios1 <- new("Scenarios")
scenarios1 <- setParamsBaseScenarios(scenarios1, horizon=5, nScenarios=10)
scenarios1 <- scenarios1 <- setRiskParamsScenarios(scenarios1, vol=.1, k=2,volStock=.2,
volRealEstate=.15, volDefault=.2, alpha=.1,beta=1, eta=.05,rho=.5, stock0=100,realEstate0=50,
liquiditySpread0=.01, defaultSpread0=.01)
}

