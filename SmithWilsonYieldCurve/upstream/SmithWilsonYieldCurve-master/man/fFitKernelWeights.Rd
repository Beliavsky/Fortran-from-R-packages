\name{fFitKernelWeights}
\alias{fFitKernelWeights}
\title{Solve for the vector xi of kernel weights}
\usage{
  fFitKernelWeights(CashflowMatrix, KernelFunctionMatrix,
    MarketValueVector, BaseZeroVector)
}
\arguments{
  \item{CashflowMatrix}{A matrix of all cashflows,
  instruments in rows, times in columns}

  \item{KernelFunctionMatrix}{A matrix of kernel function
  values}

  \item{MarketValueVector}{A vector of market values of the
  insturments}

  \item{BaseZeroVector}{A vector of "base" values for the
  zeros}
}
\description{
  Solve for the vector xi of kernel weights
}

