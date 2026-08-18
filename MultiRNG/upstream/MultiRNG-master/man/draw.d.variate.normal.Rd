\name{draw.d.variate.normal}
\alias{draw.d.variate.normal}
\title{Pseudo-Random Number Generation under Multivariate Normal Distribution}
\description{This function implements pseudo-random number generation for a multivariate normal distribution with pdf

\deqn{f(x|\mu,\Sigma)=c\exp{(-\frac{1}{2}(x-\mu)^{T}\Sigma^{-1}(x-\mu))}}

for \eqn{-\infty < x < \infty} and \eqn{c=(2\pi)^{-d/2}|\Sigma|^{-1/2}}, \eqn{\Sigma} is symmetric and positive definite, where \eqn{\mu} and \eqn{\Sigma} are the mean vector and the variance-covariance matrix, respectively.
}
\usage{
draw.d.variate.normal(no.row,d,mean.vec,cov.mat)
}
\arguments{
  \item{no.row}{Number of rows to generate.}
  \item{d}{Number of variables to generate.}
  \item{mean.vec}{Vector of means.}
  \item{cov.mat}{Variance-covariance matrix.}
}
\value{A \eqn{no.row \times d} matrix of generated data.}
\examples{
cmat<-matrix(c(1,0.2,0.3,0.2,1,0.2,0.3,0.2,1), nrow=3, ncol=3)
meanvec=c(0,3,7)
mydata=draw.d.variate.normal(no.row=1e5,d=3,mean.vec=meanvec,cov.mat=cmat)
apply(mydata,2,mean)-meanvec
cor(mydata)-cmat
}
