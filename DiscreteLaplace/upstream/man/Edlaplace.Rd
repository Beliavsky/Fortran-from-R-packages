\name{Edlaplace}
\alias{Edlaplace}
\title{
Moments of the discrete Laplace distribution
}
\description{
The function provides the expected value and the variance of the SDL, and the expectation of its absolute value.
}
\usage{
Edlaplace(p, q)
}
\arguments{
  \item{p}{
the first parameter, in \eqn{(0,1)}, of the DSL
}
  \item{q}{
the second parameter, in \eqn{(0,1)}, of the DSL
}
}
\details{
\eqn{E(X;p,q)=\frac{1}{1-p}-\frac{1}{1-q}=\frac{p}{1-p}-\frac{q}{1-q}},

\eqn{E(|X|;p,q)=\frac{q(1-p)^2+p(1-q)^2}{(1-qp)(1-q)(1-p)}},

\eqn{V(X;p,q)=\frac{1}{(1-p)^2(1-q)^2}[\frac{q(1-p)^3(1+q)+p(1-q)^3(1+p)}{1-pq}-(p-q)^2]}
}
\value{
A list of three items:
  \item{E1}{expected value}
  \item{E1a}{expectation of the absolute value}
  \item{V}{variance}
}
\references{
T. J. Kozubowski, S. Inusah (2006) A skew Laplace distribution on integers, \emph{Annals of the Institute of Statistical Mathematics}, 58: 555-571
}
\author{
Alessandro Barbiero, Riccardo Inchingolo
}
\seealso{
\code{\link{ddlaplace}}
}
\examples{
# ex.1
p<-0.5
q<-0.4
Edlaplace(p, q)
# ex.2
p<-0.1
q<-0.9
Edlaplace(p, q)
}
\keyword{distribution}
