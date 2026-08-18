\name{draw.multivariate.laplace}
\alias{draw.multivariate.laplace}
\title{Pseudo-Random Number Generation under Multivariate Laplace Distribution}
\description{This function implements pseudo-random number generation for a multivariate Laplace (double exponential) distribution with pdf

\deqn{f(x|\mu,\Sigma,\gamma)=c\exp(-((x-\mu)^{T}\Sigma^{-1}(x-\mu))^{\gamma/2})}

for \eqn{-\infty < x < \infty} and \eqn{c=\frac{\gamma\Gamma(d/2)}{2\pi^{d/2}\Gamma(d/\gamma)}|\Sigma|^{-1/2}}, \eqn{\Sigma} is symmetric and positive definite, where \eqn{\mu}, \eqn{\Sigma}, and \eqn{\gamma} are the mean vector, the variance-covariance matrix, and the shape parameter, respectively.
}
\usage{
draw.multivariate.laplace(no.row,d,gamma,mu,Sigma)
}
\arguments{
  \item{no.row}{Number of rows to generate.}
  \item{d}{Number of variables to generate.}
  \item{gamma}{Shape parameter.}
  \item{mu}{Vector of means.}
  \item{Sigma}{Variance-covariance matrix.}
}
\value{A \eqn{no.row \times d} matrix of generated data.}
\seealso{\code{\link{generate.point.in.sphere}}}
\references{
Ernst, M. D. (1998). A multivariate generalized Laplace distribution. \emph{Computational Statistics}, \bold{13}, 227-232.
}
\examples{
cmat<-matrix(c(1,0.2,0.3,0.2,1,0.2,0.3,0.2,1), nrow=3, ncol=3)
mu.vec=c(0,3,7)
mydata=draw.multivariate.laplace(no.row=1e5,d=3,gamma=2,mu=mu.vec,Sigma=cmat)
\donttest{
apply(mydata,2,mean)-mu.vec
cor(mydata)-cmat
}
}


