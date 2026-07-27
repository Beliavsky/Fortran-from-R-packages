\name{Joint.Wright}
\alias{Joint.Wright}
\title{ A Joint Version of Wight's Rank and Sign Test }
\description{
  This function returns joint or multiple version of Wright's rank and sign tests.
  The test takes the maximum value of the individual rank or sign tests, in the same manner as Chow-Denning test
}
\usage{
Joint.Wright(y, kvec)
}
\arguments{
  \item{y}{ a vector of time series, typically  financial return }
  \item{kvec}{ a vector of holding periods }
}

\value{
  \item{Holding.Period }{holding periods used}
  \item{JR1 }{Joint test based on R1 statistics}
  \item{JR2 }{Joint test based on R2 statistics}
  \item{JS1 }{Joint test based on S1 statistics}
}
\references{ 
Belaire-Franch G, Contreras D. Ranks and signs-based multiple variance ratio tests, Working paper, University of Valencia 2004. 

Kim, J. H. and Shamsuddin, A., 2008, Are Asian Stock Markets Efficient? Evidence from New Multiple Variance Ratio Tests, Journal of Empirical Fiance 15(8), 518-532.
}
\author{Jae H. Kim}
\examples{
data(exrates)
y <- exrates$ca                           
nob <- length(y)
r <- log(y[2:nob])-log(y[1:(nob-1)])    
kvec <- c(2,5,10)
Joint.Wright(r,kvec) 
}
\keyword{  htest  }