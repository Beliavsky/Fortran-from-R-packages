\name{rStock}
\alias{rStock}
\title{rStock}
\usage{
  rStock(horizon, nScenarios, ZC, vol, k, volStock, stock0,
    rho)
}
\arguments{
  \item{horizon}{Horizon of projection}

  \item{nScenarios}{Number of scenarios}

  \item{ZC}{ZC rate input}

  \item{vol}{Volatility for short rates}

  \item{k}{k for rates in vasicek model}

  \item{volStock}{Volatility}

  \item{stock0}{Initial value for stock}

  \item{rho}{Correlation between stock and short rates}
}
\description{
  Direct stock generation. Object creation is managed
  internally.
}
\examples{
data(ZC)
rStock(horizon=10, nScenarios=7, ZC=ZC, vol=.1, k=2, volStock=.2, stock0=100, rho=.5)
}

