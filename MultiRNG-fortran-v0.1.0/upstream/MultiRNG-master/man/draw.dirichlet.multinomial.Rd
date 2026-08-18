\name{draw.dirichlet.multinomial}
\alias{draw.dirichlet.multinomial}
\title{Pseudo-Random Number Generation under Dirichlet-Multinomial Distribution}
\description{This function implements pseudo-random number generation for a Dirichlet-multinomial distribution. This is a mixture distribution that is multinomial with parameter \eqn{\theta} that is a realization of a random variable having a Dirichlet distribution with shape vector \eqn{\alpha}. \eqn{N} is the sample size and \eqn{\beta} is a common scale parameter.
}
\usage{
draw.dirichlet.multinomial(no.row,d,alpha,beta,N)
}
\arguments{
  \item{no.row}{Number of rows to generate.}
  \item{d}{Number of variables to generate.}
  \item{alpha}{Vector of shape parameters.}
  \item{beta}{Scale parameter common to \eqn{d} variables.}
  \item{N}{Sample size.}
}
\value{A \eqn{no.row \times d} matrix of generated data.}
\seealso{\code{\link{draw.dirichlet}}, \code{\link{draw.multinomial}}}
\examples{
alpha.vec=c(1,3,4,4) ; N=3
mydata=draw.dirichlet.multinomial(no.row=1e5,d=4,alpha=alpha.vec,beta=2, N=3)
apply(mydata,2,mean)-N*alpha.vec/sum(alpha.vec)
}
