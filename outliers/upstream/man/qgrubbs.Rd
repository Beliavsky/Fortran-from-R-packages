\name{qgrubbs}
\alias{qgrubbs}
\alias{pgrubbs}
\title{ Calculate critical values and p-values for Grubbs tests }
\description{

This function is designed to calculate critical values for Grubbs tests
for outliers detecting and to approximate p-values reversively.

}
\usage{
qgrubbs(p, n, type = 10, rev = FALSE)
pgrubbs(q, n, type = 10)
}
\arguments{
  \item{p}{ vector of probabilities. }
  \item{q}{ vector of quantiles. }
  \item{n}{ sample size. }
  \item{type}{ Integer value indicating test variant. 10 is a test for one outlier (side is
	detected automatically and can be reversed by \code{opposite} parameter). 11 is a test
	for two outliers on opposite tails, 20 is test for two outliers in one tail. }
  \item{rev}{ if set to TRUE, function \code{qgrubbs} acts as \code{pgrubbs}. }
}
\details{

The critical values for test for one outlier is calculated according to approximations given by
Pearson and Sekar (1936). The formula is simply reversed to obtain p-value.

The values for two outliers test (on opposite sides) are calculated according to David, Hartley,
and Pearson (1954). Their formula cannot be rearranged to obtain p-value, thus such values are
obtained by \code{\link{uniroot}}.

For test checking presence of two outliers at one tail, the tabularized distribution (Grubbs, 1950)
is used, and approximations of p-values are interpolated using \code{\link{qtable}}.

}
\value{

A vector of quantiles or p-values.

}
\references{ 

Grubbs, F.E. (1950). Sample Criteria for testing outlying observations. Ann. Math. Stat. 21, 1, 27-58.

Pearson, E.S., Sekar, C.C. (1936). The efficiency of statistical tools and a criterion for the rejection of outlying observations. Biometrika, 28, 3, 308-320.

David, H.A, Hartley, H.O., Pearson, E.S. (1954). The distribution of the ratio, in a single normal sample, of range to standard deviation. Biometrika, 41, 3, 482-493.
}
\author{ Lukasz Komsta }
\seealso{ \code{\link{grubbs.test}}}
\keyword{distribution}
