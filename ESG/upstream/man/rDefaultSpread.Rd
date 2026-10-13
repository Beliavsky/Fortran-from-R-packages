\name{rDefaultSpread}
\alias{rDefaultSpread}
\title{rDefaultSpread}
\usage{
  rDefaultSpread(horizon, nScenarios, defaultSpread0,
    volDefault, alpha, beta)
}
\arguments{
  \item{horizon}{Horizon of projection}

  \item{nScenarios}{Number of scenarios}

  \item{defaultSpread0}{Initial value for default spread}

  \item{volDefault}{Volatility}

  \item{alpha}{alpha for LMN model}

  \item{beta}{beta Volatility for LMN model}
}
\description{
  Direct default spread generation. Object creation is
  managed internally.
}
\examples{
rDefaultSpread(horizon=5, nScenarios=8, defaultSpread0=.01, volDefault=.2, alpha=.1, beta=1)
}

