\name{rpois}
\alias{rpois.ccomp}
\alias{rmultinom.ccomp}
\title{Simulate count compositions without overdispersion}
\description{
  Generates multinomial or multi-Poission random variates based on an
  Aitchison composition.
}
\usage{
rpois.ccomp(n,p,lambda)
rmultinom.ccomp(n,p,N)
          }
\arguments{
  \item{n}{number of datasets to be simulated}
  \item{p}{The composition representing the probabilites/portions of the
  individual parts}
  \item{lambda}{scalar or vector giving the expected total count}
  \item{N}{scalar or vector giving the total count}
}
\value{
  a random dataset ccount dataset
}
\details{
  A count composition is a realisation of a multinomial or multivariate
  Poisson distribution. 
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{rnorm.ccomp}}
}
\examples{
p <- acomp(c(3,3,3))
rpois.ccomp(10,p,40)
rmultinom.ccomp(10,p,40)

}
\keyword{distribution}
