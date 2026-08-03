\name{ArchTest}
\alias{ArchTest}

\title{
  ARCH LM Test 
}
\description{
  Lagrange Multiplier (LM) test for autoregressive conditional
  heteroscedasticity (ARCH) 
}
\usage{
ArchTest (x, lags=12, demean = FALSE) 
}
\arguments{
  \item{x}{numeric vector} 
  \item{lags}{positive integer number of lags} 
  \item{demean}{
    logical:  If TRUE, remove the mean before computing the test
    statistic.  
  } 
}
\details{
  Computes the Lagrange multiplier test for conditional
  heteroscedasticity of Engle (1982), as described by Tsay (2005,
  pp. 101-102).

  This is provided for compatibility with 'archTest' in the S-Plus
  script in Tsay (p. 102).  
}
\value{
  an object of class 'htest'  
}

\seealso{
  \code{\link{AutocorTest}}
}
\author{Bernhard Pfaff}
\examples{
data(m.intc7303)
intcLM <- ArchTest(log(1+as.numeric(m.intc7303)), lag=12)
# Matches answer on Tsay (p. 102) 
}
\keyword{ts}
