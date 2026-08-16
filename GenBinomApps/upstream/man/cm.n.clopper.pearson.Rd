\name{cm.n.clopper.pearson}
\alias{cm.n.clopper.pearson}
\title{
Required Sample Size - Countermeasure Model
}
\description{
Provides the required sample size with respect to the extended upper Clopper-Pearson limit for a failure model, where countermeasures are introduced. }
\usage{
cm.n.clopper.pearson(p, size, cm.effect, alpha = 0.1, uniroot.lower = k + 1,
  uniroot.upper = 1e+100, uniroot.tol = 1e-10, uniroot.maxiter = 1e+05)
}
\arguments{
  \item{p}{
target probability.}
  \item{size}{
vector of the number of failures for each type.
}
  \item{cm.effect}{
vector of the success probabilities to solve a failure for each type. Corresponds to the probabilities \eqn{p_{i}}{pi} of a generalized binomial distribution.
}
  \item{alpha}{
significance level for the \eqn{(1-\alpha)\cdot 100}{(1-alpha)* 100}\%  confidence level (default \eqn{\alpha=0.1}{alpha=0.1}).
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
Provides the required sample size with respect to the extended upper Clopper-Pearson limit. It applies for the case that failures are tackled by countermeasures. That means countermeasures with different effectivities for each failure type are introduced. See the references for further information.
}
\value{
The value for the required sample size.
}
\references{
D.Kurz, H.Lewitschnig, J.Pilz, \emph{Decision-Theoretical Model for Failures which are Tackled by Countermeasures}, IEEE Transactions on Reliability, Vol. 63, No. 2, June 2014.
}

\seealso{
 \code{\link[stats]{uniroot}},\code{\link{dgbinom}},\code{\link{cm.clopper.pearson.ci}},\code{\link{n.clopper.pearson}}
}
\examples{
## target failure probability p=0.00001, 2 failures divided in 2 failure types k1=1, k2=1.
## 2 countermeasures with effectivities p1=0.5, p2=0.8
cm.n.clopper.pearson(0.00001,size=c(1,1),cm.effect=c(0.5,0.8))
# 365299
}