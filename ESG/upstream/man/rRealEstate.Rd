\name{rRealEstate}
\alias{rRealEstate}
\title{rRealEstate}
\usage{
  rRealEstate(horizon, nScenarios, ZC, vol, k,
    volRealEstate, realEstate0)
}
\arguments{
  \item{horizon}{Horizon of projection}

  \item{nScenarios}{Number of scenarios}

  \item{ZC}{ZC rate input}

  \item{vol}{Volatility for short rates}

  \item{k}{k for rates in vasicek model}

  \item{volRealEstate}{Volatility}

  \item{realEstate0}{Initial value for real estate}
}
\description{
  Direct real estate generation. Object creation is managed
  internally.
}
\examples{
data(ZC)
rRealEstate(horizon=5, nScenarios=10, ZC=ZC, vol=.1, k=2, volRealEstate=.15, realEstate0=50)
}

