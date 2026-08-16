\encoding{UTF-8}

\name{rargus}
\alias{rargus}
\alias{dargus}
\alias{pargus}
\alias{qargus}

\title{ 
  The Argus distribution
}

\description{
Density, distribution function, quantile function and random generation for the Argus distribution with parameter \code{chi}.
\code{rargus()} uses very efficient random variate generation methods. Inversion or the Ratio of Uniforms method can be selected. 
}

\usage{
rargus(n=length(chi) , chi, method = c("inversion","RoU") )
dargus(x, chi, log = FALSE)
pargus(x, chi, lower = TRUE, log.p = FALSE)
qargus(p, chi, lower = TRUE, log.p = FALSE)}

\arguments{
  \item{n}{Number of observations}
  \item{chi}{parameter vector \code{chi>0}}
  \item{method}{random variate generation method to be used (details are below)}
  \item{x}{x-value for the density or CDF}
  \item{log}{If \code{TRUE} the logarithm of the density is returned.}
  \item{lower}{If \code{FALSE} 1-CDF or its inverse is returned}
  \item{log.p}{\code{log(p)} the logarithm of the probability is used or returned}
  \item{p}{probability of the quantile}
}

\details{
  The Argus distribution has the density (pdf) which is proportional to
  \deqn{f(x) = x \sqrt{1-x^2}\exp(-0.5\chi^2(1-x^2)),  \mbox{for} 0 \leq x \leq 1 \mbox{and} \chi > 0}{%
        f(x) = x \sqrt(1-x^2)\exp(-0.5\chi^2(1-x^2)),  for 0 <= x <= 1 and \chi > 0}


  The generators used for \code{rargus()} are very efficient also for the varying parameter case,
  i.e. when the vector \code{chi} has length \code{n}.
  
  \code{method = "inversion"}: transforms the uniform variate into an argus variate
   using a close to exact approximation of the inverse CDF. It is also well suited for the 
   varying parameter case as the theorem that an argus variate can be represented as a transformed 
   truncated Gamma(1.5) variate allows to express the inverse CDF of the argus distribution as a 
   simple transform of the inverse CDF of the Gamma(1.5) distribution.
   The use of the \code{pinv.new()} function of the \code{Runuran} package makes that evaluation very fast. 
   When loading the argus-package \code{pinv.new()} calculates and stores the required tables.
   
   \code{method = "RoU"}: uses the Ratio of Uniforms method which requires 2 or more uniform variates
   to generate one argus variate. It is also a bit slower than the inversion method and is added here
   mainly for the case that the \code{Runuran}-package is not available.
}

\value{
\code{rargus} creates a random sample of size \code{n}. 
 
\code{dargus} gives the density, \code{dargus} give the CDF and \code{qargus} gives the quantile function.
The length of the result for these 3 functions is the maximum of the lengths of the numerical arguments that
are recycled to the length of the result. Only the first elements of the logical arguments are used.
}

\author{
  Wolfgang \enc{Hörmann}{Hoermann} and
  Christoph Baumgarten
}

\references{%
  Christoph Baumgarten: Random Variate Generation by Fast Numerical Inversion in the Varying Parameter Case. 
}

%% \seealso{}

\examples{
library(argus)
# pdf, cdf and quantile function
dargus((0:10)/10,.6)
qargus(0.9,2)
pargus(qargus(0.9,2),2)

y<-rargus(n=1.e5,chi=0.3,method="inversion")
hist(y,breaks=100,main="Argus pdf with chi=0.3",freq=FALSE)
lines(xv<-seq(0,1,0.001),dargus(xv,chi=0.3),lwd=0.2,col=2)

n=1.e5
chiv <- 1. +runif(n)*5 # chi-values in (1,5)
system.time({set.seed(123);y<-rargus(chi=chiv)})
set.seed(123);u<- runif(n)
mean(abs(pargus(y,chiv)-u))

# random variate generation with the inversion method
y<-rargus(n=1.e5,chi=0.3,method="inversion")
hist(y,breaks=100,main="Argus pdf with chi=0.3",freq=FALSE)
lines(xv<-seq(0,1,0.001),dargus(xv,chi=0.3),lwd=0.2,col=2)
# using Ratio of Uniforms
y<-rargus(n=1.e5,chi=3,method="RoU")
hist(y,breaks=100,main="Argus pdf with chi=3",freq=FALSE)
lines(xv<-seq(0,1,0.001),dargus(xv,chi=3),lwd=2,col=2)
# generating for different chi values
y<-rargus(n=100,chi=runif(100)*2,method="RoU")
}

\keyword{package}
\keyword{distribution}
\keyword{datagen}
