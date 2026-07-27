\name{Wright.crit}
\alias{Wright.crit}
\title{ Critical Values for Wright's rank and sign tests }
\description{
  This function returns critical values of Wright's tests based on the simulation method detailed in Wright (2000)
}
\usage{
Wright.crit(n, k, nit)
}
\arguments{
  \item{n}{ sample size }
  \item{k}{ holding period, a scalar }
  \item{nit}{ number of iterations }
}

\value{
  \item{Holding.Period }{holding period used}
  \item{R1.crit }{Critical values for the R1 statistic}
   \item{R2.crit }{Critical values for the R2 statistic}
   \item{S1.crit }{Critical values for the S1 statistic}
}
\references{ WRIGHT,J.H.,2000,Alternative Variance-Ratio Tests Using Ranks and Signs, Journal of Business & Economic Statistics, 18, 1-9. }
\author{ Jae H. Kim }

\examples{
Wright.crit(n=10,k=2,nit=50)  

# nit is set to 50 for fast execution in the example. 
# nit=10000 is recommended as in Wright (2000)
}
\keyword{  htest  }