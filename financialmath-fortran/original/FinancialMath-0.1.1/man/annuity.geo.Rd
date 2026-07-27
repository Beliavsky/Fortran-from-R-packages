\arguments{
  \item{pv}{present value of the annuity}
  \item{fv}{future value of the annuity}
  \item{n}{number of payments/periods for the annuity}
  \item{p}{amount of the first payment}
  \item{k}{payment growth rate per period}
  \item{i}{nominal interest rate convertible ic times per year}
  \item{ic}{interest conversion frequency per year}
  \item{pf}{the payment frequency- number of payments/periods per year}
  \item{imm}{option for annuity immediate or annuity due, default is immediate (TRUE)}
  \item{plot}{option to display a time diagram of the payments}
}
\description{Solves for the present value, future value, number of payments/periods, amount of the first payment, the payment growth rate, and/or the interest rate for a geometrically growing annuity.  It can also plot a time diagram of the payments.}
\name{annuity.geo}
\title{Geometric Annuity}
\alias{annuity.geo}
\keyword{annuity}
\keyword{geometric}
\note{At least one of pv, fv, n, pmt, k, or i must be NA (unknown).

pv and fv cannot both be specified, at least one must be NA (unknown).}
\value{Returns a matrix of the input variables and calculated unknown variables.}
\examples{
annuity.geo(pv=NA,fv=100,n=10,p=9,k=.02,i=NA,ic=2,pf=.5,plot=TRUE)

annuity.geo(pv=NA,fv=128,n=5,p=NA,k=.04,i=.03,pf=2)}
\usage{annuity.geo(pv=NA,fv=NA,n=NA,p=NA,k=NA,i=NA,ic=1,pf=1,imm=TRUE,plot=FALSE)}
\seealso{
  \code{\link{annuity.arith}}

  \code{\link{annuity.level}}

  \code{\link{perpetuity.arith}}

  \code{\link{perpetuity.geo}}

  \code{\link{perpetuity.level}}
}

\details{
Effective Rate of Interest: \eqn{eff.i=(1+\frac{i}{ic})^{ic}-1}

\eqn{j=(1+eff.i)^{\frac{1}{pf}}-1}

\eqn{fv=pv*(1+j)^n}

Annuity Immediate:

j != k: \eqn{pv=p*\frac{1-(\frac{1+k}{1+j})^n}{j-k}}

j = k: \eqn{pv=p*\frac{n}{1+j}}

Annuity Due:

j != k: \eqn{pv=p*\frac{1-(\frac{1+k}{1+j})^n}{j-k}*(1+j)}

j = k: \eqn{pv=p*n}
}
