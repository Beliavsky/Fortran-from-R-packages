\arguments{
  \item{pv}{present value}
  \item{pmt}{value of level payments}
  \item{i}{nominal interest rate convertible ic times per year}
  \item{ic}{interest conversion frequency per year}
  \item{pf}{the payment frequency- number of payments per year}
  \item{imm}{option for perpetuity immediate or annuity due, default is immediate (TRUE)}
}
\name{perpetuity.level}
\alias{perpetuity.level}
\title{Level Perpetuity}
\description{Solves for the present value, interest rate, or the amount of the payments for a level perpetuity.}
\usage{perpetuity.level(pv=NA,pmt=NA,i=NA,ic=1,pf=1,imm=TRUE)}
\author{Kameron Penn and Jack Schmidt}
\keyword{perpetuity}
\keyword{level}
\note{One of pv, pmt, or i must be NA (unknown).}
\value{Returns a matrix of the input variables and calculated unknown variables.}
\examples{
perpetuity.level(pv=100,pmt=NA,i=.05,ic=1,pf=2,imm=TRUE)

perpetuity.level(pv=100,pmt=NA,i=.05,ic=1,pf=2,imm=FALSE)}
\seealso{
  \code{\link{perpetuity.arith}}

  \code{\link{perpetuity.geo}}

  \code{\link{annuity.arith}}

  \code{\link{annuity.geo}}

  \code{\link{annuity.level}}
}

\details{
Effective Rate of Interest: \eqn{eff.i=(1+\frac{i}{ic})^{ic}-1}

\eqn{j=(1+eff.i)^{\frac{1}{pf}}-1}

Perpetuity Immediate:

\eqn{pv=pmt*{a_{\left. {\overline {\, \infty \,}}\! \right |j}}=\frac{pmt}{j}}

Perpetuity Due:

\eqn{pv=pmt*{\ddot {a}_{\left. {\overline {\, \infty \,}}\! \right |j}}=\frac{pmt}{j}*(1+i)}

}
