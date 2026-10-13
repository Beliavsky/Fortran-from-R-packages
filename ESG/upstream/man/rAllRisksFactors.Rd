\name{rAllRisksFactors}
\alias{rAllRisksFactors}
\title{rAllRisksFactors}
\usage{
  rAllRisksFactors(horizon, nScenarios, ZC, vol, k,
    volStock, stock0, rho, volRealEstate, realEstate0, eta,
    liquiditySpread0, defaultSpread0, volDefault, alpha,
    beta)
}
\arguments{
  \item{horizon}{Horizon of projection}

  \item{nScenarios}{Number of scenarios}

  \item{ZC}{ZC rate input}

  \item{vol}{Volatility for short rates}

  \item{k}{k for rates in vasicek model}

  \item{volStock}{Volatility for stock}

  \item{stock0}{Initial value for stock}

  \item{rho}{Correlation between stock and short rates}

  \item{volRealEstate}{Volatility for real estate}

  \item{realEstate0}{Initial value for real estate}

  \item{eta}{eta Volatility for LMN model}

  \item{liquiditySpread0}{Initial value for liquidity
  spread}

  \item{defaultSpread0}{Initial value for default spread}

  \item{volDefault}{Volatilty for default spread}

  \item{alpha}{alpha for LMN model}

  \item{beta}{beta Volatility for LMN model}
}
\description{
  Direct generation for all risk factors. Object creation
  is managed internally.
}
\examples{
data(ZC)
rAllRisksFactors(horizon=5, nScenarios=10, ZC, vol=.1, k=2, volStock=.2, stock0=100, rho=.5,
volRealEstate=.15, realEstate0=50, eta=.05, liquiditySpread0=.01, defaultSpread0=.01,
volDefault=.2, alpha=.1, beta=1)
}

