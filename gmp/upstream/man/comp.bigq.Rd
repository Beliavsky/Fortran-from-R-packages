\name{Bigq}

\alias{==.bigq}
\alias{!=.bigq}
\alias{<.bigq}
\alias{<=.bigq}
\alias{>.bigq}
\alias{>=.bigq}
\alias{sign.bigq}
\title{Relational Operators}
\description{
  Binary operators which allow the comparison of values in atomic vectors.
}
\usage{
\method{sign}{bigq}(x)

\method{<}{bigq}(e1, e2)
\method{<=}{bigq}(e1, e2)
\method{==}{bigq}(e1, e2)
\method{>=}{bigq}(e1, e2)
\method{>}{bigq}(e1, e2)
\method{!=}{bigq}(e1, e2)
}
\arguments{
  \item{x, e1, e2}{Object or vector of class \code{\link{bigq}}}
}

\examples{
  x <- as.bigq(8000,21)
  x  < 2 * x
  }

\keyword{arith}

