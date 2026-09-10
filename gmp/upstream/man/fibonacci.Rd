\name{lucnum}
\alias{fibnum}
\alias{fibnum2}
\alias{lucnum}
\alias{lucnum2}

\title{Compute Fibonacci and Lucas numbers}
\description{
  fibnum  compute n-th Fibonacci number.
  fibnum2 compute (n-1)-th and n-th Fibonacci number.
  lucnum  compute n-th lucas number.
  lucnum2 compute (n-1)-th and n-th lucas number.

  Fibonacci numbers are define by: \eqn{F_n=F_{n-1}+F_{n-2}}{Fn=Fn-1 + Fn-2}
  Lucas numbers are define by: \eqn{L_n=F_n+2F_{n-1}}{Ln=Fn + 2Fn-1}

}
\usage{
fibnum(n)
fibnum2(n)
lucnum(n)
lucnum2(n)
}

\arguments{
  \item{n}{Integer}
}

\value{
Fibonacci numbers and Lucas number.
}

\references{The GNU MP Library, see \url{https://gmplib.org}}
\author{Antoine Lucas}



\examples{
fibnum(10)
fibnum2(10)
lucnum(10)
lucnum2(10)

}
\keyword{arith}

