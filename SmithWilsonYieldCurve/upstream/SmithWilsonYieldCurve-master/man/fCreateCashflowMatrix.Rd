\name{fCreateCashflowMatrix}
\alias{fCreateCashflowMatrix}
\title{Returns the matrix of cashflows for the list of instruments}
\usage{
  fCreateCashflowMatrix(dfInstruments)
}
\arguments{
  \item{dfInstruments}{A set of market instruments as a
  dataframe with columns Type, Tenor, Frequency and Rate
  with Type in (LIBOR, SWAP, BOND), Tenor the instrument maturity
  in years and rate the rate per annum}
}
\description{
  Returns the matrix of cashflows for the list of
  instruments
}

