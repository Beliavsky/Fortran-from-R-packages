\name{draw.d.variate.t}
\alias{draw.d.variate.t}
\title{Pseudo-Random Number Generation under Multivariate t Distribution}
\description{This function implements pseudo-random number generation for a multivariate \eqn{t} distribution with pdf

\deqn{f(x|\mu, \Sigma, \nu)=c\left(1+\frac{1}{\nu}(x-\mu)^{T}\Sigma^{-1}(x-\mu)\right)^{-(\nu+d)/2}}

for \eqn{-\infty < x < \infty} and \eqn{c=\frac{\Gamma((\nu+d)/2)}{\Gamma(\nu/2)(\nu\pi)^{d/2}}|\Sigma|^{-1/2}}, \eqn{\Sigma} is symmetric and positive definite, \eqn{\nu>0}, where \eqn{\mu}, \eqn{\Sigma}, and \eqn{\nu} are the mean vector, the variance-covariance matrix, and the degrees of freedom, respectively.

}
\usage{
draw.d.variate.t(dof,no.row,d,mean.vec,cov.mat)
}
\arguments{
  \item{dof}{Degrees of freedom.}
  \item{no.row}{Number of rows to generate.}
  \item{d}{Number of variables to generate.}
  \item{mean.vec}{Vector of means.}
  \item{cov.mat}{Variance-covariance matrix.}
}
\value{A \eqn{no.row \times d} matrix of generated data.}
\examples{
cmat<-matrix(c(1,0.2,0.3,0.2,1,0.2,0.3,0.2,1), nrow=3, ncol=3)
meanvec=c(0,3,7)
mydata=draw.d.variate.t(dof=5,no.row=1e5,d=3,mean.vec=meanvec,cov.mat=cmat)
apply(mydata,2,mean)-meanvec
cor(mydata)-cmat
}
