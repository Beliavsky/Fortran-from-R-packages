\name{payoff}
\alias{payoff}
\title{
Multi-Asset Option Payoff Calculator
}
\description{
\code{payoff} calculates the per-unit option payoff for digital, best-of, and worst-of multi-asset options.
}
\usage{
payoff(payType, pcFlag, strike, S)
}
\arguments{
  \item{payType}{
    case; if 0, digital payoff, if 1, best-of payoff, if 2, worst-of payoff.
  }
  \item{pcFlag}{
    case vector; if 0, call, if 1, put.
  }
  \item{strike}{
    vector; option strikes.
  }
  \item{S}{
    list containing the vectors of spatial grid points associated with each underlying. Vector sizes of underlying spatial grid points need not be equal.
  }
}
\value{
\code{payoff} returns an array of the unit option values at each point spanned by the list of underlying vectors. Dimension of array is inhereted from \code{S}.
}
\author{
Michael Eichenberger and Carlo Rosa
}
\examples{
# payoff of a dual-asset digital call with strikes at 100 and 90.
S <- list(seq(0, 500, by=1), seq(0, 500, by=1))
payoff(0, c(0, 0), c(100, 90), S)
}
