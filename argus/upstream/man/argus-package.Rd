\encoding{UTF-8}
\name{argus-package}
\alias{argus-package}
\alias{argus}
\docType{package}

\title{
  Generator and density for the Argus distribution
}

\description{
  Random variate generation, density, CDF and quantile function 
  for the Argus distribution. Especially, it includes for random variate generation a
  flexible inversion method that is also fast in the varying parameter case.
  A Ratio-of-Uniforms method is provided as second alternative.
}

\details{
  The Argus distribution has the density (pdf) which is proportional to
  \deqn{f(x) = x \sqrt{1-x^2}\exp(-0.5\chi^2(1-x^2)),  \mbox{for} 0 \leq x \leq 1 \mbox{and} \chi > 0}{%
        f(x) = x \sqrt(1-x^2)\exp(-0.5\chi^2(1-x^2)),  for 0 <= x <= 1 and \chi > 0}

  \tabular{ll}{
    Package: \tab argus          \cr
    License: \tab GPL 2 or later\cr
  }

  Package \pkg{argus} provides the routines:
  \describe{
    \item{\code{\link{rargus}}}{generates argus distributed random variates.}
    \item{\code{\link{dargus}}}{computes the density of the argus distribution.}
    \item{\code{\link{pargus}}}{computes the CDF of the argus distribution.}
    \item{\code{\link{qargus}}}{computes the quantile function (ie. the inverse CDF )of the argus distribution.}
  }

}

\author{
  Wolfgang \enc{Hörmann}{Hoermann} & Christoph Baumgarten
}

\references{%
Christoph Baumgarten: Random Variate Generation by Fast Numerical Inversion in the Varying Parameter Case. 
}

\examples{
## Evaluate the pdf (density)
dargus(c(0.1,0.5,0.9), chi=0.3)
## Evaluate the CDF for different chi values
pargus(c(0.1,0.5,0.9), chi=c(0.3,1.3,2.3)) 
## Evaluate the quantile function
pargus(c(0.1,0.5,0.9), chi=4.5)
## Draw a random sample
rargus(n=10, chi=0.3)
## compare histogram and density
system.time(y<-rargus(1.e5,chi=2.5))
hist(y,breaks=100,freq=FALSE)
lines(x<-seq(0,1,1.e-3),dargus(x,2.5),col=2,lwd=2)
}

\keyword{package}
\keyword{distribution}
\keyword{datagen}
