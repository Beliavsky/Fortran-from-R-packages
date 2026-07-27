\name{JWright.crit}
\alias{JWright.crit}
\title{Critical Values for the joint versions of Wright's rank and sign tests  }
\description{
  This function runs a simulation to calculate the critical values of the joint versions of Wright's tests.
}
\usage{
JWright.crit(n, kvec, nit)
}
\arguments{
  \item{n}{ sample size }
  \item{kvec}{ holding period vector }
  \item{nit}{ number of iterations }
}

\value{
  \item{Holding.Period }{holding period used}
  \item{JR1.crit }{Critical values for the joint R1 statistic}
   \item{JR2.crit }{Critical values for the joint R2 statistic}
   \item{JS1.crit }{Critical values for the joint S1 statistic}
}
\references{

Belaire-Franch G, Contreras D. Ranks and signs-based multiple variance ratio tests, Working paper, University of Valencia 2004. 

Kim, J. H. and Shamsuddin, A., 2008, Are Asian Stock Markets Efficient? Evidence from New Multiple Variance Ratio Tests, Journal of Empirical Fiance 15(8), 518-532.
}

\author{ Jae H. Kim }

\examples{

kvec <- c(2,5,10)
JWright.crit(n=100,kvec,nit=50)

# nit is set to 50 for fast execution in the example. 
# nit=10000 is recommended as in Wright (2000)
}
\keyword{  htest  }