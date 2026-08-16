\name{Edlaplace2}
\alias{Edlaplace2}
\title{
First- and second-order moments of ADSL
}
\description{
First- and second-order moment of the ADSL distribution.
}
\usage{
Edlaplace2(p, q)
}
\arguments{
  \item{p}{
the first parameter \eqn{p}, in \eqn{(0,1)}, of the ADSL
}
  \item{q}{
the first parameter \eqn{q}, in \eqn{(0,1)}, of the ADSL
}
}
\details{
For the ADSL distribution,

\eqn{E(X;p,q)=\frac{\log q}{\log(pq)}\frac{p}{1-p}-\frac{\log p}{\log(pq)}\frac{1}{1-q}}

and

\eqn{E(X^2;p,q)=\frac{\log q}{\log(pq)}\frac{p(1+p)}{(1-p)^2}+\frac{\log p}{\log(pq)}\frac{1+q}{(1-q)^2}}
}
\value{
A list containing the first- and the second-order moments of the ADSL distribution, \code{E1} and \code{E2}.
}
\references{
A. Barbiero, An alternative discrete Laplace distribution, \emph{Statistical Methodology}, 16: 47-67
}
\author{
Alessandro Barbiero, Riccardo Inchingolo
}
\seealso{
\code{\link{estdlaplace2}}, \code{\link{loss}}, \code{\link{ddlaplace2}}
}
\examples{
Edlaplace2(p=0.3, q=0.3)
Edlaplace2(p=0.3, q=0.6)
Edlaplace2(p=0.6, q=0.3)
Edlaplace2(p=0.6, q=0.6)
}
\keyword{distribution}
