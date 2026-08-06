\name{fCreateKernelMatrix}
\alias{fCreateKernelMatrix}
\title{Create the matrix of kernel functions}
\usage{
  fCreateKernelMatrix(times, fKernel)
}
\arguments{
  \item{times}{a vector of cashflow times}

  \item{fKernel}{a kernel to apply (a function of times x
  times returning a matrix )}
}
\description{
  Creates a J x J matrix [ w(u_i,u_j) ] where J is the
  number of cashflow times in the calibration set
}

