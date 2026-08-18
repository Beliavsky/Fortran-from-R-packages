\name{draw.wishart}
\alias{draw.wishart}
\title{Pseudo-Random Number Generation under Wishart Distribution}
\description{This function implements pseudo-random number generation for a Wishart distribution with pdf

\deqn{f(x|\nu,\Sigma)=(2^{\nu d/2}\pi^{d(d-1)/4}\prod_{i=1}^{d}\Gamma((\nu+1-i)/2))^{-1}|\Sigma|^{-\nu/2}|x|^{(\nu-d-1)/2}\exp(-\frac{1}{2}tr(\Sigma^{-1}x))}

\eqn{x} is positive definite, \eqn{\nu \geq d}, and \eqn{\Sigma} is symmetric and positive definite, where \eqn{\nu} and \eqn{\Sigma} are positive definite and the scale matrix, respectively.
}
\usage{
draw.wishart(no.row,d,nu,sigma)
}
\arguments{
  \item{no.row}{Number of rows to generate.}
  \item{d}{Number of variables to generate.}
  \item{nu}{Degrees of freedom.}
  \item{sigma}{Scale matrix.}
}
\value{A \eqn{no.row \times d^2} matrix of Wishart deviates in the form of rows.To obtain the Wishart matrix, convert each row to a matrix where rows are filled first.}
\seealso{\code{\link{draw.d.variate.normal}}}
\examples{
\donttest{
mymat<-matrix(c(1,0.2,0.3,0.2,1,0.2,0.3,0.2,1), nrow=3, ncol=3)
draw.wishart(no.row=1e5,d=3,nu=5,sigma=mymat)
}
}

