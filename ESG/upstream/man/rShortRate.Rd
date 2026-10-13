\name{rShortRate}
\alias{rShortRate}
\title{rShortRate}
\usage{
  rShortRate(horizon, nScenarios, ZC, vol, k)
}
\arguments{
  \item{horizon}{Horizon of projection}

  \item{nScenarios}{Number of scenarios}

  \item{ZC}{ZC rate input}

  \item{vol}{Volatility for short rates}

  \item{k}{k for rates in vasicek model}
}
\description{
  Direct short rate generation. Object creation is managed
  internally.
}
\examples{
data(ZC)
rShortRate(horizon=15, nScenarios=10, ZC=ZC, vol=.1, k=2)
}

