\name{amort.table}
\alias{amort.table}
\title{Amortization Table}
\description{Produces an amortization table for paying off a loan while also solving for either the number of payments, loan amount, or the payment amount.  In the amortization table the payment amount, interest paid, principal paid, and balance of the loan are given for each period.  If n ends up not being a whole number, outputs for the balloon payment, drop payment and last regular payment are provided. The total interest paid, and total amount paid is also given.  It can also plot the percentage of each payment toward interest vs. period.}
\usage{amort.table(Loan=NA,n=NA,pmt=NA,i,ic=1,pf=1,plot=FALSE)}
\author{Kameron Penn and Jack Schmidt}
\arguments{
  \item{Loan}{loan amount}
  \item{n}{the number of payments/periods}
  \item{pmt}{value of level payments}
  \item{i}{nominal interest rate convertible ic times per year}
  \item{ic}{interest conversion frequency per year}
  \item{pf}{the payment frequency- number of payments per year}
  \item{plot}{tells whether or not to plot the percentage of each payment toward interest vs. period}
}
\examples{
amort.table(Loan=1000,n=2,i=.005,ic=1,pf=1)

amort.table(Loan=100,pmt=40,i=.02,ic=2,pf=2,plot=FALSE)

amort.table(Loan=NA,pmt=102.77,n=10,i=.005,plot=TRUE)}
\note{Assumes that payments are made at the end of each period.

One of n, Loan, or pmt must be NA (unknown).

If pmt is less than the amount of interest accumulated in the first period, then the function will stop because the loan will never be paid off due to the payments being too small.

If pmt is greater than the loan amount plus interest accumulated in the first period, then the function will stop because one payment will pay off the loan.}
\keyword{amortization}
\value{A list of two components.
\item{Schedule}{A data frame of the amortization schedule.}
\item{Other}{A matrix of the input variables and other calculated variables.}
}
\seealso{
  \code{\link{amort.period}}

  \code{\link{annuity.level}}
}
\details{
Effective Rate of Interest: \eqn{eff.i=(1+\frac{i}{ic})^{ic}-1}

\eqn{j=(1+eff.i)^{\frac{1}{pf}}-1}

\eqn{Loan=pmt*{a_{\left. {\overline {\, n \,}}\! \right |j}}}

Balance at the end of period t: \eqn{B_t=pmt*{a_{\left. {\overline {\, n-t \,}}\! \right |j}}}

Interest paid at the end of period t: \eqn{i_t=B_{t-1}*j}

Principal paid at the end of period t: \eqn{p_t=pmt-i_t}

Total Paid\eqn{=pmt*n}

Total Interest Paid\eqn{=pmt*n-Loan}

If \eqn{n=n^*+k} where \eqn{n^*} is an integer and \eqn{0<k<1}:

Last regular payment (at period \eqn{n^*}) \eqn{=pmt*{s_{\left. {\overline {\, k \,}}\! \right |j}}}

Drop payment (at period \eqn{n^*+1}) \eqn{=Loan*(1+j)^{n^*+1}-pmt*{s_{\left. {\overline {\, n^* \,}}\! \right |j}}}

Balloon payment (at period \eqn{n^*}) \eqn{=Loan*(1+j)^{n^*}-pmt*{s_{\left. {\overline {\, n^* \,}}\! \right |j}}+pmt}
}
