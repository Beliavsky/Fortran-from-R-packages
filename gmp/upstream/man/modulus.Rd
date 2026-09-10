\name{modulus}
\alias{modulus}
\alias{modulus<-}
\alias{modulus.bigz}
\alias{modulus<-.bigz}
\title{Modulus of a Big Integer}

\description{
  The modulus of a \code{\link{bigz}} number \eqn{a} is
  \dQuote{unset} when \eqn{a} is a regular integer, \eqn{a \in Z}).
  Or the modulus can be set to \eqn{m} which means
  \eqn{a \in Z/\,m\cdot Z}{a in Z/mZ}), i.e., all arithmetic with
  \eqn{a} is performed \sQuote{modulo m}.
}
\usage{
modulus(a)
modulus(a) <- value
}
\arguments{
  \item{a}{\R object of class \code{"\link{bigz}"}}
  \item{value}{integer number or object of class \code{"\link{bigz}"}.}
}

\examples{
x <- as.bigz(24)
modulus(x) # NULL, i.e. none

# x element of Z/31Z :
modulus(x) <- 31
x+x  #  48 |-> (17 \%\% 31)
10*x # 240 |-> (23 \%\% 31)
x31 <- x

# reset modulus to "none":
modulus(x) <- NA; x; x. <- x
x <- x31
modulus(x) <- NULL; x

stopifnot(identical(x,            as.bigz(24)), identical(x, x.),
          identical(modulus(x31), as.bigz(31)))
}
\keyword{arith}
