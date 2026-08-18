\name{draw.multinomial}
\alias{draw.multinomial}
\title{Pseudo-Random Number Generation under Multivariate Multinomial Distribution}
\description{This function implements pseudo-random number generation for a multivariate multinomial distribution with pdf

\deqn{f(x|\theta_{1},...,\theta_{d})=\frac{N!}{\prod x_{j}!}\prod_{j=1}^{d}\theta_{j}^{x_{j}}}

for \eqn{0<\theta_{j}<1}, \eqn{x_{j}\geq 0}, and \eqn{\sum_{j=1}^{d}x_{j}=N}, where \eqn{\theta_{1},...,\theta_{d}} are cell probabilities and \eqn{N} is the size.
}
\usage{
draw.multinomial(no.row,d,theta,N)
}
\arguments{
  \item{no.row}{Number of rows to generate.}
  \item{d}{Number of variables to generate.}
  \item{theta}{Vector of cell probabilities.}
  \item{N}{Sample Size. Must be at least 2.}
}
\value{A \eqn{no.row \times d} matrix of generated data.}
\examples{
theta.vec=c(0.3,0.3,0.25,0.15) ; N=4
mydata=draw.multinomial(no.row=1e5,d=4,theta=c(0.3,0.3,0.25,0.15),N=4)
apply(mydata,2,mean)-N*theta.vec
}
