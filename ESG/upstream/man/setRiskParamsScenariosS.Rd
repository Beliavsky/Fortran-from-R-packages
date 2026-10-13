\name{setRiskParamsScenariosS}
\alias{setRiskParamsScenariosS}
\alias{setRiskParamsScenariosS,Scenarios-method}
\title{setRiskParamsScenariosS method}
\arguments{
  \item{vol}{Volatility for rates in vasicek model}

  \item{k}{k for rates in vasicek model}

  \item{volStock}{Volatility for UL in Black & Scholes
  model}

  \item{stock0}{UL initial value}

  \item{rho}{Correlation between stock and short rates}
}
\description{
  Set risk parameters related to short rates in a Scenarios
  object (these parameters are contained in a
  [\code{\linkS4class{ParamsScenarios}}] sub-object)
}

