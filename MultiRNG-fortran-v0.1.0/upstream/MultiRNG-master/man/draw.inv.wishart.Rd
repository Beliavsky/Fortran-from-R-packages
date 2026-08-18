\name{draw.inv.wishart}
\alias{draw.inv.wishart}
\title{Pseudo-Random Number Generation under Inverted Wishart Distribution}
\description{This function implements pseudo-random number generation for an inverted Wishart distribution with pdf

\deqn{f(x|\nu,\Sigma)=(2^{\nu d/2}\pi^{d(d-1)/4}\prod_{i=1}^{d}\Gamma((\nu+1-i)/2))^{-1}|\Sigma|^{\nu/2}|x|^{-(\nu+d+1)/2}\exp(-\frac{1}{2}tr(\Sigma x^{-1}))}

\eqn{x} is positive definite, \eqn{\nu \geq d}, and \eqn{\Sigma^{-1}} is symmetric and positive definite, where \eqn{\nu} and \eqn{\Sigma^{-1}} are the degrees of freedom and the inverse scale matrix, respectively.
}
\usage{
draw.inv.wishart(no.row,d,nu,inv.sigma)
}
\arguments{
  \item{no.row}{Number of rows to generate.}
  \item{d}{Number of variables to generate.}
  \item{nu}{Degrees of freedom.}
  \item{inv.sigma}{Inverse scale matrix.}
}
\value{A \eqn{no.row \times d^2} matrix ofcontaining Wishart deviates in the form of rows. To obtain the Inverted-Wishart matrix, convert each row to a matrix where rows are filled first.}
\seealso{\code{\link{draw.wishart}}}
\examples{
\donttest{
mymat<-matrix(c(1,0.2,0.3,0.2,1,0.2,0.3,0.2,1), nrow=3, ncol=3)
draw.inv.wishart(no.row=1e5,d=3,nu=5,inv.sigma=mymat)
}
}

