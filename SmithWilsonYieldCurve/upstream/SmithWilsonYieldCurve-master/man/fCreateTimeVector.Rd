\name{fCreateTimeVector}
\alias{fCreateTimeVector}
\title{Extract a vector of cashflow times in years from a list of instruments}
\usage{
  fCreateTimeVector(dfInstruments)
}
\arguments{
  \item{dfInstruments}{A dataframe of instuments with at
  least columns Type and Tenor}
}
\description{
  Assumes that LIBOR tenor is in days, with 365 days per
  year. Assumes that SWAPs are semi-annual Returns a vector
  of all unique cashflow times in years
}

