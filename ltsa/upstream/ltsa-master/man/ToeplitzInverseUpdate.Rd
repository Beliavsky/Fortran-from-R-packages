\name{ToeplitzInverseUpdate}
\alias{ToeplitzInverseUpdate}
\title{ Inverse of Toeplitz matrix of order n+1 given inverse of order n }
\description{
Let G be a Toeplitz matrix of order n and with (i,j)-element, r[Abs[i-j]].
So the first row of G may be written (r[0],...,r[n-1]).  Suppose the next
element in the sequence is r[n].  Then the inverse of the Toeplitz matrix
whose first row is (r[0],...,r[n]) may be obtained either using
ToeplitzInverseUpdate or directly using TrenchInverse.  ToeplitzInverseUpdate
is somewhat faster. 
}
\usage{
ToeplitzInverseUpdate(GI, r, rnew)
}
\arguments{
  \item{GI}{ inverse of Toeplitz matrix G of order n }
  \item{r}{ first row of G , ie r[0],...,r[n-1]}
  \item{rnew}{ next element, r[n] }
}
\details{
Although this update requires \eqn{O(n^2)}{O(n^2)} flops, the same as TrenchInverse,
it is somewhat faster in practice.
}
\value{
inverse matrix of order n+1
}
\references{
Graybill, F.A. (1983). Matrices with Applications in Statistics.

McLeod, A.I., Yu, Hao, Krougly, Zinovi L.  (2007).
Algorithms for Linear Time Series Analysis,
Journal of Statistical Software.

}
\author{ A.I. McLeod }
\seealso{ 
  \code{\link{TrenchInverse}}
}
\examples{
#In this example we compute the update inverse directly and using ToeplitzInverseUpdate and
#compare the result.
phi<-0.8
sde<-30
n<-30
r<-arima.sim(n=30,list(ar=phi),sd=sde)
r<-phi^(0:(n-1))/(1-phi^2)*sde^2
n1<-25
G<-toeplitz(r[1:n1])
GI<-solve(G) #could also use TrenchInverse
GIupdate<-ToeplitzInverseUpdate(GI,r[1:n1],r[n1+1])
GIdirect<-solve(toeplitz(r[1:(n1+1)]))
ERR<-sum(abs(GIupdate-GIdirect))
ERR
}
\keyword{ ts }
\keyword{array }

