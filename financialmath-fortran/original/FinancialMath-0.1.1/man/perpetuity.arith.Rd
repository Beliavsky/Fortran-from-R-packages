\name{perpetuity.arith}
\title{Arithmetic Perpetuity}
\alias{perpetuity.arith}
\description{Solves for the present value, amount of the first payment, the payment increment amount per period, or the interest rate for an arithmetically growing perpetuity.}
\usage{perpetuity.arith(pv=NA,p=NA,q=NA,i=NA,ic=1,pf=1,imm=TRUE)}
\author{Kameron Penn and Jack Schmidt}
\arguments{
  \item{pv}{present value of the annuity}
  \item{p}{amount of the first payment}
  \item{q}{payment increment amount per period}
  \item{i}{nominal interest rate convertible ic times per year}
  \item{ic}{interest conversion frequency per year}
  \item{pf}{the payment frequency- number of payments per year}
  \item{imm}{option for annuity immediate or annuity due, default is immediate (TRUE)}
}
\examples{
perpetuity.arith(100,p=1,q=.5,i=NA,ic=1,pf=1,imm=TRUE)

perpetuity.arith(pv=NA,p=1,q=.5,i=.07,ic=1,pf=1,imm=TRUE)

perpetuity.arith(pv=100,p=NA,q=1,i=.05,ic=.5,pf=1,imm=FALSE)}
\note{One of pv, p, q, or i must be NA (unknown).}
\keyword{perpetuity}
\keyword{arithmetic}
\value{Returns a matrix of input variables, and calculated unknown variables.}
\seealso{
  \code{\link{perpetuity.geo}}

  \code{\link{perpetuity.level}}

  \code{\link{annuity.arith}}

  \code{\link{annuity.geo}}

  \code{\link{annuity.level}}
}

\details{
Effective Rate of Interest: \eqn{eff.i=(1+\frac{i}{ic})^{ic}-1}

\eqn{j=(1+eff.i)^{\frac{1}{pf}}-1}

Perpetuity Immediate:

\eqn{pv=\frac{p}{j}+\frac{q}{j^2}}

Perpetuity Due:

\eqn{pv=(\frac{p}{j}+\frac{q}{j^2})*(1+j)}
}
