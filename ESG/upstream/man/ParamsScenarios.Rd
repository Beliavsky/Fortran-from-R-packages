\name{ParamsScenarios}
\alias{ParamsScenarios}
\alias{ParamsScenarios-class}
\title{ParamsScenarios class}
\description{
  This class is a container for all the risk related
  parameters. It is used as a parameter for the Scenarios
  class.
}
\details{
  \describe{ \item{horizon}{Horizon for the projection (in
  years)}

  \item{nScenarios}{Number of scenarios}

  \item{vol}{Volatility for rates in vasicek model}

  \item{k}{k for rates in vasicek model}

  \item{volStock}{Volatility for UL in Black & Scholes
  model}

  \item{volRealEstate}{Volatility for real estate in Black
  & Scholes model}

  \item{stock0}{Stock initial value}

  \item{realEstate0}{Real estate initial value}

  \item{volDefault}{Volatility for LMN model}

  \item{alpha}{alpha for LMN model}

  \item{beta}{beta for LMN model}

  \item{eta}{eta for LMN model}

  \item{liquiditySpread0}{Initial liquidity for LMN model}

  \item{defaultSpread0}{Initial default spread for LMN
  model}

  \item{rho}{Correlation between stock and short rates} }
}

