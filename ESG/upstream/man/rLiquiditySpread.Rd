\name{rLiquiditySpread}
\alias{rLiquiditySpread}
\title{rLiquiditySpread}
\usage{
  rLiquiditySpread(horizon, nScenarios, eta,
    liquiditySpread0)
}
\arguments{
  \item{horizon}{Horizon of projection}

  \item{nScenarios}{Number of scenarios}

  \item{eta}{eta Volatility for LMN model}

  \item{liquiditySpread0}{Initial value for liquidity
  spread}
}
\description{
  Direct liquidity spread generation. Object creation is
  managed internally.
}
\examples{
rLiquiditySpread(horizon=5, nScenarios=15, eta=.05, liquiditySpread0=.01)
}

