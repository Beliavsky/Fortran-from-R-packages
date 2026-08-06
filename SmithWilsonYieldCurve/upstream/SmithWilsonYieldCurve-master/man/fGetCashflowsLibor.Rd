\name{fGetCashflowsLibor}
\alias{fGetCashflowsLibor}
\title{Gets the cashflow schedule for a LIBOR agreement}
\usage{
  fGetCashflowsLibor(dfInstrument)
}
\arguments{
  \item{dfInstrument}{A set of market instruments as a
  dataframe with columns Type, Tenor and Rate with Type in
  (LIBOR, SWAP), Tenor the instrument maturity in years and
  rate the rate per annum}
}
\description{
  Gets the cashflow schedule for a LIBOR agreement
}

