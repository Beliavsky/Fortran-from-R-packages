\name{iFI}
\alias{iFI}
\title{
Inverse of Fisher Information matrix
}
\description{
Inverse of Fisher Information matrix for the DSL.
}
\usage{
iFI(p, q)
}
\arguments{
  \item{p}{
first parameter \eqn{p}, in \eqn{(0,1)}, of the DSL
}
  \item{q}{
second parameter \eqn{q}, in \eqn{(0,1)}, of the DSL
}
}
\value{
The inverse of Fisher Information matrix. Take a look at the references for more details.
}
\references{
T. J. Kozubowski, S. Inusah (2006) A skew Laplace distribution on integers, \emph{Annals of the Institute of Statistical Mathematics}, 58: 555-571
}
\author{
Alessandro Barbiero, Riccardo Inchingolo
}
\seealso{
\code{\link{iFI2}}
}
\examples{
p <- 0.2
q <- 0.8
iFI(p, q)
}
\keyword{distribution,htest}
