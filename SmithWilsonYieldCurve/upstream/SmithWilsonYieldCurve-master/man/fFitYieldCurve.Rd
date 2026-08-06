\name{fFitYieldCurve}
\alias{fFitYieldCurve}
\title{Constructs the ZCB function based on the given market inputs and a specific kernel and base function}
\usage{
  fFitYieldCurve(TimesVector, CashflowMatrix,
    MarketValueVector, fKernel, fBase)
}
\arguments{
  \item{TimesVector}{A vector of all cashflow times}

  \item{CashflowMatrix}{A matrix of all cashflows,
  instruments in rows, times in columns}

  \item{MarketValueVector}{A vector of market values of the
  insturments}

  \item{fKernel}{a function of two times used as the Kernel
  "basis" function}

  \item{fBase}{a function giving the base level of the
  curve}
}
\value{
  a list comprising elements: a function of time which
  gives the ZCB price to that time
}
\description{
  Constructs the ZCB function based on the given market
  inputs and a specific kernel and base function
}

