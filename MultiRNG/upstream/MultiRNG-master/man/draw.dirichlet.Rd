\name{draw.dirichlet}
\alias{draw.dirichlet}
\title{Pseudo-Random Number Generation under Multivariate Beta (Dirichlet) Distribution}
\description{This function implements pseudo-random number generation for a multivariate beta (Dirichlet) distribution with pdf

\deqn{f(x|\alpha_{1},...,\alpha_{d})=\frac{\Gamma(\sum_{j=1}^{d}\alpha_{j})}{\prod_{j=1}^{d}\Gamma(\alpha_{j})} \prod_{j=1}^{d}x_{j}^{\alpha_{j}-1}}

for \eqn{\alpha_{j}>0}, \eqn{x_{j}\geq 0}, and \eqn{\sum_{j=1}^{d}x_{j}=1}, where \eqn{\alpha_{1},...,\alpha_{d}} are the shape parameters and \eqn{\beta} is a common scale paramter.
}
\usage{
draw.dirichlet(no.row,d,alpha,beta)
}
\arguments{
  \item{no.row}{Number of rows to generate.}
  \item{d}{Number of variables to generate.}
  \item{alpha}{Vector of shape parameters.}
  \item{beta}{Scale parameter common to \eqn{d} variables.}
}
\value{A \eqn{no.row \times d} matrix of generated data.}
\examples{
alpha.vec=c(1,3,4,4)
mydata=draw.dirichlet(no.row=1e5,d=4,alpha=alpha.vec,beta=2)
apply(mydata,2,mean)-alpha.vec/sum(alpha.vec)
}
