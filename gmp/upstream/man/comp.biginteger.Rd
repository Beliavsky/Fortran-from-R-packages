\name{Relational_Operator}

\alias{==.bigz}
\alias{!=.bigz}
\alias{<.bigz}
\alias{<=.bigz}
\alias{>.bigz}
\alias{>=.bigz}
\alias{sign.bigz}
\title{Relational Operators}
\description{
  Binary operators which allow the comparison of values in atomic vectors.
}
\usage{
\method{sign}{bigz}(x)
\method{==}{bigz}(e1, e2)
\method{<}{bigz} (e1, e2)
\method{>=}{bigz}(e1, e2)
}
\arguments{
  \item{x, e1, e2}{\R object (vector or matrix-like) of class \code{"\link{bigz}"}.}
}
\seealso{\code{\link{mod.bigz}} for arithmetic operators.
}
\examples{
  x <- as.bigz(8000)
  x ^ 300 < 2 ^x

  sign(as.bigz(-3:3))
  sign(as.bigq(-2:2, 7))
}

\keyword{arith}

