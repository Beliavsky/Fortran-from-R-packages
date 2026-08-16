\name{estdlaplace}
\alias{estdlaplace}
\title{
Sample estimation for the DSL
}
\description{
The function provides the maximum likelihood estimates for the parameters of the DSL and the estimate of the inverse of the Fisher information matrix. The method of moments estimates of \eqn{p} and \eqn{q} coincide with the maximum likelihood estimates.
}
\usage{
estdlaplace(x)
}
\arguments{
  \item{x}{
a vector of observations from the DSL
}
}
\details{
See the reference.
If
\eqn{\bar{x}^{+}=\frac{1}{n}\sum_{i=1}^n x_i^{+}}, \eqn{\bar{x}^{-}=\frac{1}{n}\sum_{i=1}^n x_i^{-}}
where \eqn{x^{+}} and \eqn{x^{-}} are the positive and the negative parts of \eqn{x}, respectively: \eqn{x^{+}=x} if \eqn{x\geq 0} and zero otherwise, \eqn{x^{-}=(-x)^{+}}, then

\eqn{\hat{q}=\frac{2\bar{x}^{-}(1+\bar{x})}{1+2\bar{x}^{-}\bar{x}+\sqrt{1+4\bar{x}^{-}\bar{x}^{+}}}}, 
\eqn{\hat{p}=\frac{\hat{q}+\bar{x}(1-\hat{q})}{1+\bar{x}(1-\hat{q})}}

when \eqn{\bar{x}\geq 0} and

\eqn{\hat{p}=\frac{2\bar{x}^{+}(1-\bar{x})}{1-2\bar{x}^{+}\bar{x}+\sqrt{1+4\bar{x}^{-}\bar{x}^{+}}}}, 
\eqn{\hat{q}=\frac{\hat{p}-\bar{x}(1-\hat{p})}{1-\bar{x}(1-\hat{p})}}

when \eqn{\bar{x}\leq 0}.
}
\value{
A list comprising
\item{hatp}{estimate of \eqn{p}}
\item{hatq}{estimate of \eqn{q}}
\item{hatSigma}{estimate of the inverse of the Fisher information matrix}
}
\references{
T. J. Kozubowski, S. Inusah (2006) A skew Laplace distribution on integers, \emph{Annals of the Institute of Statistical Mathematics}, 58: 555-571
}
\author{
Alessandro Barbiero, Riccardo Inchingolo
}
\seealso{
\code{\link{ddlaplace}}
}
\examples{
p<-0.6
q<-0.3
n<-20
x<-rdlaplace(n, p, q)
est<-estdlaplace(x)
est[1]
est[2]
est[3]
# increase n
n<-100
x<-rdlaplace(n, p, q)
est<-estdlaplace(x)
est[1]
est[2]
est[3]
# swap the parameters
x<-rdlaplace(n, q, p)
est<-estdlaplace(x)
est[1]
est[2]
est[3]
}
\keyword{distribution}
\keyword{htest}

