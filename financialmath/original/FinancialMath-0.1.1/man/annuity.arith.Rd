\name{annuity.arith}
\title{Arithmetic Annuity}
\alias{annuity.arith}
\description{Solves for the present value, future value, number of payments/periods, amount of the first payment, the payment increment amount per period, and/or the interest rate for an arithmetically growing annuity.  It can also plot a time diagram of the payments.}
\usage{annuity.arith(pv=NA,fv=NA,n=NA,p=NA,q=NA,i=NA,ic=1,pf=1,imm=TRUE,plot=FALSE)}
\author{Kameron Penn and Jack Schmidt}
\arguments{
  \item{pv}{present value of the annuity}
  \item{fv}{future value of the annuity}
  \item{n}{number of payments/periods}
  \item{p}{amount of the first payment}
  \item{q}{payment increment amount per period}
  \item{i}{nominal interest frequency convertible ic times per year}
  \item{ic}{interest conversion frequency per year}
  \item{pf}{the payment frequency- number of payments per year}
  \item{imm}{option for annuity immediate or annuity due, default is immediate (TRUE)}
  \item{plot}{option to display a time diagram of the payments}
}
\examples{
annuity.arith(pv=NA,fv=NA,n=20,p=100,q=4,i=.03,ic=1,pf=2,imm=TRUE)

annuity.arith(pv=NA,fv=3000,n=20,p=100,q=NA,i=.05,ic=3,pf=2,imm=FALSE)
}
\note{At least one of pv, fv, n, p, q, or i must be NA (unknown).

pv and fv cannot both be specified, at least one must be NA (unknown).}
\keyword{annuity}
\keyword{arithmetic}
\value{Returns a matrix of the input variables, and calculated unknown variables.}
\seealso{
  \code{\link{annuity.geo}}

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

\eqn{pv=p*{a_{\left. {\overline {\, n \,}}\! \right |j}}+q* \frac{{a_{\left. {\overline {\, n \,}}\! \right |j}}-n*(1+j)^{-n}}{j}}

Annuity Due:

\eqn{pv=(p*{a_{\left. {\overline {\, n \,}}\! \right |j}}+q* \frac{{a_{\left. {\overline {\, n \,}}\! \right |j}}-n*(1+j)^{-n}}{j})*(1+i)}

}
