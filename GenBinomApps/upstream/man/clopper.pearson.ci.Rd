\name{clopper.pearson.ci}
\alias{clopper.pearson.ci}
\title{Clopper-Pearson Confidence Interval}
\description{ Computing upper, lower or two-sided Clopper-Pearson confidence limits for a given confidence  level.
}
\usage{
clopper.pearson.ci(k, n, alpha = 0.1, CI = "upper")
}
\arguments{
  \item{k}{
number of failures/successes.
}
  \item{n}{
number of trials.
}
  \item{alpha}{
significance level for the \eqn{(1-\alpha)\cdot 100}{(1-alpha)* 100}\%  confidence level (default \eqn{\alpha=0.1}{alpha=0.1}).
}
  \item{CI}{
indicates the kind of the confidence interval, options: "upper" (default), "lower",  "two.sided". 
}
}
\details{
Computes the confidence limits for the \eqn{p} of a binomial distribution.
Confidence intervals are obtained by the definition of Clopper and Pearson.
The two-sided interval for \eqn{k=0} is \eqn{(0,1-(\alpha/2)^{1/n})}{(0,1-(alpha/2)^(1/n))}, for \eqn{k=n} it is \eqn{((\alpha/2)^{1/n},1)}{((alpha/2)^(1/n),1)}.
}
\value{
A data frame containing the kind of the confidence interval, upper and lower limits and the used significance level \code{alpha}.
}
\references{
D.Kurz, H.Lewitschnig, J.Pilz, \emph{Decision-Theoretical Model for Failures which are Tackled by Countermeasures}, IEEE Transactions on Reliability, Vol. 63, No. 2, June 2014.

Thulin, Mans, \emph{The cost of using exact confidence intervals for a binomial proportion}, Electronic Journal of Statistics, vol. 8, pp. 817-840, 2014.

C.J.Clopper and E.S. Pearson, \emph{The use of confidence or fiducial limits illustrated in the case of the binomial}, Biometrika, vol. 26, pp. 404-413, 1934.
}

\examples{
clopper.pearson.ci(5,100000,alpha=0.05)
# Confidence.Interval = upper
# Lower.limit = 0
# Upper.limit = 0.0001051275
# alpha = 0.05

clopper.pearson.ci(5,100000,CI="two.sided")
# Confidence.Interval =  two.sided
# Lower.limit = 1.97017e-05
# Upper.limit = 0.0001051275
# alpha = 0.1
}
