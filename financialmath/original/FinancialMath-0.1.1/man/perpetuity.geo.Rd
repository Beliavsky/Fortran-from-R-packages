\arguments{
  \item{pv}{present value}
  \item{p}{amount of the first payment}
  \item{k}{payment growth rate per period}
  \item{i}{nominal interest rate convertible ic times per year}
  \item{ic}{interest conversion frequency per year}
  \item{pf}{the payment frequency- number of payments and periods per year}
  \item{imm}{option for perpetuity immediate or due, default is immediate (TRUE)}
}
\description{Solves for the present value, amount of the first payment, the payment growth rate, or the interest rate for a geometrically growing perpetuity. }
\name{perpetuity.geo}
\alias{perpetuity.geo}
\title{Geometric Perpetuity}
\keyword{perpetuity}
\keyword{geometric}
\note{One of pv, p, k, or i must be NA (unknown).}
\value{Returns a matrix of the input variables and calculated unknown variables.}
\examples{
perpetuity.geo(pv=NA,p=5,k=.03,i=.04,ic=1,pf=1,imm=TRUE)

perpetuity.geo(pv=1000,p=5,k=NA,i=.04,ic=1,pf=1,imm=FALSE)}
\usage{perpetuity.geo(pv=NA,p=NA,k=NA,i=NA,ic=1,pf=1,imm=TRUE)}
\seealso{
  \code{\link{perpetuity.arith}}

  \code{\link{perpetuity.level}}

  \code{\link{annuity.arith}}

  \code{\link{annuity.geo}}

  \code{\link{annuity.level}}
}

\details{
Effective Rate of Interest: \eqn{eff.i=(1+\frac{i}{ic})^{ic}-1}

\eqn{j=(1+eff.i)^{\frac{1}{pf}}-1}

Perpetuity Immediate:

j > k: \eqn{pv=\frac{p}{j-k}}

Perpetuity Due:

j > k: \eqn{pv=\frac{p}{j-k}*(1+j)}

}
