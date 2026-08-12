\name{dirichletrnd}
\alias{dirichletrnd}
\title{Dirichlet generator}
\usage{
  dirichletrnd(a, n) 
}
\description{
  Random generation for the Dirichlet distribution
}

\arguments{
\item{a}{numeric vector for the concentration parameters}
\item{n}{number of observations}
}

\value{dirichletrnd generates \code{n} random observations from a Dirichlet distribution}
	

\author{Tim Benham, Qibin Duan, Dirk P. Kroese, Benoit Liquet}

\references{Kroese D.P., Taimre T., Botev Z.I. (2011), \emph{Handbook of Monte Carlo Methods}, John Wiley & Sons.}

\examples{
## Generation from the Dirichlet distribution 
## with parameter a=(1,2,3,4,5)

set.seed(12345)
a <- 1:5
n <- 10

y <- dirichletrnd(a,n)
y
}


