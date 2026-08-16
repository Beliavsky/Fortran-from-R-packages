\name{n.clopper.pearson}
\alias{n.clopper.pearson}
\title{
Required Sample Size 
}
\description{
Provides the required sample size with respect to the one-sided upper Clopper-Pearson limit. }

\usage{
n.clopper.pearson(k, p, alpha = 0.1, uniroot.lower = k + 1, uniroot.upper = 1e+100, 
uniroot.maxiter = 1e+05, uniroot.tol = 1e-10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{
number of failures.
}
  \item{p}{
target probability.
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
Provides the required sample size with respect to the upper Clopper-Pearson limit for a given target failure probability at a certain confidence level.
}
\value{
The value for the required sample size.
}
\references{
D.Kurz, H.Lewitschnig, J.Pilz, \emph{Decision-Theoretical Model for Failures which are Tackled by Countermeasures}, IEEE Transactions on Reliability, Vol. 63, No. 2, June 2014.
}




\seealso{
 \code{\link[stats]{uniroot}},\code{\link{clopper.pearson.ci}}
}
\examples{
## target failure probability p=0.0002, 8 failures
n.clopper.pearson(8,0.0002)
# 64972
}

