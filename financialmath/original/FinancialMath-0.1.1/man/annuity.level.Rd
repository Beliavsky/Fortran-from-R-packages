\arguments{
  \item{pv}{present value of the annuity}
  \item{fv}{future value of the annuity}
  \item{n}{number of payments/periods}
  \item{pmt}{value of the level payments}
  \item{i}{nominal interest rate convertible ic times per year}
  \item{ic}{interest conversion frequency per year}
  \item{pf}{the payment frequency- number of payments/periods per year}
  \item{imm}{option for annuity immediate or annuity due, default is immediate (TRUE)}
  \item{plot}{option to display a time diagram of the payments}
}
\description{Solves for the present value, future value, number of payments/periods, interest rate, and/or the amount of the payments for a level annuity.  It can also plot a time diagram of the payments.}
\name{annuity.level}
\alias{annuity.level}
\title{Level Annuity}
\keyword{level}
\keyword{annuity}
\note{At least one of pv, fv, n, pmt, or i must be NA (unknown).

pv and fv cannot both be specified, at least one must be NA (unknown).}
\value{Returns a matrix of the input variables and calculated unknown variables.}
\examples{
annuity.level(pv=NA,fv=101.85,n=10,pmt=8,i=NA,ic=1,pf=1,imm=TRUE)

annuity.level(pv=80,fv=NA,n=15,pf=2,pmt=NA,i=.01,imm=FALSE)}
\usage{
annuity.level(pv=NA,fv=NA,n=NA,pmt=NA,i=NA,ic=1,pf=1,imm=TRUE,plot=FALSE)}
\seealso{
  \code{\link{annuity.arith}}

  \code{\link{annuity.geo}}

  \code{\link{perpetuity.arith}}

  \code{\link{perpetuity.geo}}

  \code{\link{perpetuity.level}}
}
\details{
Effective Rate of Interest: \eqn{eff.i=(1+\frac{i}{ic})^{ic}-1}

\eqn{j=(1+eff.i)^{\frac{1}{pf}}-1}

Annuity Immediate:

\eqn{pv=pmt*{a_{\left. {\overline {\, n \,}}\! \right |j}}=pmt*\frac{1-(1+j)^{-n}}{j}}

\eqn{fv=pmt*{s_{\left. {\overline {\, n \,}}\! \right |j}}=pmt*{a_{\left. {\overline {\, n \,}}\! \right |j}}*(1+j)^n}

Annuity Due:

\eqn{pv=pmt*{\ddot {a}_{\left. {\overline {\, n \,}}\! \right |j}}=pmt*{a_{\left. {\overline {\, n \,}}\! \right |j}}*(1+j)}

\eqn{fv=pmt*{\ddot {s}_{\left. {\overline {\, n \,}}\! \right |j}}=pmt*{a_{\left. {\overline {\, n \,}}\! \right |j}}*(1+j)^{n+1}}

}
