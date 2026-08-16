\name{cm.clopper.pearson.ci}
\alias{cm.clopper.pearson.ci}
\title{
Clopper-Pearson Confidence Interval for Failures Which are Tackled by Countermeasures
}
\description{
Provides the extended Clopper-Pearson confidence limits for a failure model, where countermeasures are introduced. 
}
\usage{
cm.clopper.pearson.ci(n, size, cm.effect, alpha = 0.1, CI = "upper", uniroot.lower = 0, 
 uniroot.upper = 1, uniroot.maxiter = 1e+05, uniroot.tol = 1e-10)
}
\arguments{
  \item{n}{
sample size.
}
  \item{size}{
vector of the number of  failures for each type.
}
  \item{cm.effect}{
vector of the success probabilities to solve a failure for each type. Corresponds to the probabilities \eqn{p_{i}}{pi} of a generalized binomial distribution.
}
  \item{alpha}{
significance level for the \eqn{(1-\alpha)\cdot 100}{(1-alpha)* 100}\%  confidence level (default \eqn{\alpha=0.1}{alpha=0.1}).
}
  \item{CI}{
indicates the kind of the confidence interval, options: "upper" (default), "lower",  "two.sided". 
}
  \item{uniroot.lower}{ 
The value of the \code{lower} parameter sent to \code{uniroot}. Lower bound of the interval to be searched. See \code{uniroot} for more details.
}
  \item{uniroot.upper}{
The value of the \code{upper} parameter sent to \code{uniroot}. Upper bound of the interval to be searched. See \code{uniroot} for more details.
}
  \item{uniroot.maxiter}{
The value of the \code{maxiter} parameter sent to \code{uniroot}. Maximum number of iterations. See \code{uniroot} for more details.
}
  \item{uniroot.tol}{
The value of the \code{tol} parameter sent to \code{uniroot}. Convergence tolerance. See \code{uniroot} for more details.
}
}
\details{
This is an extension of the Clopper-Pearson confidence interval, where different outcome scenarios of the random sampling are weighted by generalized binomial probabilities. The weights are the probabilities for observing \eqn{0,\dots ,k}{0,...,k} failures after the introduction of countermeasures.
Computes the confidence limits for the \eqn{p} of a binomial distribution, where \eqn{p} is the failure probability.  The failures are tackled by countermeasures for specific failure types with different effectivity. 
See the references for further information.
}
\value{
A data frame containing the kind of the confidence interval, upper and lower limits and the used significance level \code{alpha}.
}
\references{
D.Kurz, H.Lewitschnig, J.Pilz, \emph{Decision-Theoretical Model for Failures which are Tackled by Countermeasures}, IEEE Transactions on Reliability, Vol. 63, No. 2, June 2014.
}




\seealso{
 \code{\link[stats]{uniroot}}, \code{\link{dgbinom}}, \code{\link{clopper.pearson.ci}}
}
\examples{
## n=110000 tested devices, 2 failures divided in 2 failure types k1=1, k2=1.
## 2 countermeasures with effectivities p1=0.5, p2=0.8
cm.clopper.pearson.ci(110000,size=c(1,1),cm.effect=c(0.5,0.8))
# Confidence.Interval = upper
# Lower.limit = 0
# Upper.limit = 3.32087e-05
# alpha = 0.1

}

