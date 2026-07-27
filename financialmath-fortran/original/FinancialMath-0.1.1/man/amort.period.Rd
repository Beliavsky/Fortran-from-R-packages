\name{amort.period}
\alias{amort.period}
\title{Amortization Period}
\description{Solves for either the number of payments, the payment amount, or the amount of a loan. The payment amount, interest paid, principal paid, and balance of the loan are given for a specified period.}
\usage{amort.period(Loan=NA,n=NA,pmt=NA,i,ic=1,pf=1,t=1)}
\author{Kameron Penn and Jack Schmidt}
\arguments{
  \item{Loan}{loan amount}
  \item{n}{the number of payments/periods}
  \item{pmt}{value of level payments}
  \item{i}{nominal interest rate convertible ic times per year}
  \item{ic}{interest conversion frequency per year}
  \item{pf}{the payment frequency- number of payments per year}
  \item{t}{the specified period for which the payment amount, interest paid, principal paid, and loan balance are solved for}
}
\examples{
amort.period(Loan=100,n=5,i=.01,t=3)

amort.period(n=5,pmt=30,i=.01,t=3,pf=12)

amort.period(Loan=100,pmt=24,ic=1,i=.01,t=3)
}
\note{Assumes that payments are made at the end of each period.

One of n, pmt, or Loan must be NA (unknown).

If pmt is less than the amount of interest accumulated in the first period, then the function will stop because the loan will never be paid off due to the payments being too small.

If the pmt is greater than the loan amount plus interest accumulated in the first period, then the function will stop because one payment will pay off the loan.

t cannot be greater than n.}
\keyword{amortization}
\value{Returns a matrix of input variables, calculated unknown variables, and amortization figures for the given period.}
\seealso{
  \code{\link{amort.table}}
}
\details{
Effective Rate of Interest: \eqn{eff.i=(1+\frac{i}{ic})^{ic}-1}

\eqn{j=(1+eff.i)^{\frac{1}{pf}}-1}

\eqn{Loan=pmt*{a_{\left. {\overline {\, n \,}}\! \right |j}}}

Balance at the end of period t: \eqn{B_t=pmt*{a_{\left. {\overline {\, n-t \,}}\! \right |j}}}

Interest paid at the end of period t: \eqn{i_t=B_{t-1}*j}

Principal paid at the end of period t: \eqn{p_t=pmt-i_t}
}
