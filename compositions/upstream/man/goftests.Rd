\name{gof}
\alias{acompGOF.test}
\alias{acompGOF.test.list}
\alias{acompGOF.test.formula}
\alias{gsi.acompUniformityGOF.test}
\alias{acompNormalGOF.test}
\alias{acompTwoSampleGOF.test}
\title{Compositional Goodness of fit test}
\description{
  Goodness of fit tests for compositional data.
}
\usage{
acompGOF.test(x,...)
acompNormalGOF.test(x,...,method="etest")
\S3method{acompGOF.test}{formula}(formula, data,...,method="etest")
\S3method{acompGOF.test}{list}(x,...,method="etest")
gsi.acompUniformityGOF.test(x,samplesize=nrow(x)*20,R=999)
acompTwoSampleGOF.test(x,y,...,method="etest",data=NULL)
          }
\arguments{
  \item{x}{a dataset of compositions (acomp)}
  \item{y}{a dataset of compositions (acomp)}
  \item{samplesize}{number of observations in a reference sample
    specifying the distribution to compare with. Typically substantially
    larger than the sample under investigation}
  \item{R}{The number of replicates to compute the distribution
    of the test statistic}
  \item{method}{Selecting a method to be used. Currently only "etest"
    for using an energy test is supported.}
  \item{\dots}{further arguments to the methods}
  \item{formula}{an anova model formula defining groups in the dataset }
  \item{data}{unused}
}
\value{
  A classical \code{"htest"} object
  \item{data.name}{The name of the dataset as specified}
  \item{method}{a name for the test used}
  \item{alternative}{an empty string}
  \item{replicates}{a dataset of p-value distributions under the
    Null-Hypothesis got from
    nonparametric bootstrap}
  \item{p.value}{The p.value computed for this test}
}
\details{
  The compositional goodness of fit testing problem is essentially a
  multivariate goodness of fit test. However there is a lack of
  standardized multivariate goodness of fit tests in R. Some can be found in
  the \code{energy}-package.

  In principle there is only one test behind the Goodness of fit tests
  provided here, a two sample test with test statistic.
  \deqn{\frac{\sum_{ij} k(x_i,y_i)}{\sqrt{\sum_{ij} k(x_i,x_i)\sum_{ij} k(y_i,y_i)}}}
  The idea behind that statistic is to measure the cos of an angle
  between the distributions in a scalar product given by
  \deqn{
  (X,Y)=E[k(X,Y)]=E[\int K(x-X)K(x-Y) dx]
  }
  where k and K are Gaussian kernels with different spread. The bandwith
  is actually the standarddeviation of k.\cr
  The other goodness of fit tests against a specific distribution are
  based on estimating the parameters of the distribution, simulating a
  large dataset of that distribution and apply the two sample goodness
  of fit test. 
  
  For the moment, this function covers: two-sample tests, uniformity tests and additive
  logistic normality tests. Dirichlet distribution tests will be included soon.
}
\section{Missing Policy}{
  Up to now the tests can not handle missings. 
  }

\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}

\seealso{
  \code{\link{fitDirichlet}},\code{\link{rDirichlet}}, \code{\link{runif.acomp}},
  \code{\link{rnorm.acomp}}, 
}

\references{
Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
    Data} Monographs on Statistics and Applied Probability. Chapman &
  Hall Ltd., London (UK). 416p.\cr
}
\examples{
\dontrun{
x <- runif.acomp(100,4)
y <- runif.acomp(100,4)

erg <- acompTwoSampleGOF.test(x,y)
#continue
erg
unclass(erg)
erg <- acompGOF.test(x,y)


x <- runif.acomp(100,4)
y <- runif.acomp(100,4)
dd <- replicate(1000,acompGOF.test(runif.acomp(100,4),runif.acomp(100,4))$p.value)
hist(dd)

dd <- replicate(1000,acompGOF.test(runif.acomp(20,4),runif.acomp(100,4))$p.value)
hist(dd)
dd <- replicate(1000,acompGOF.test(runif.acomp(10,4),runif.acomp(100,4))$p.value)

hist(dd)
dd <- replicate(1000,acompGOF.test(runif.acomp(10,4),runif.acomp(400,4))$p.value)
hist(dd)
dd <- replicate(1000,acompGOF.test(runif.acomp(400,4),runif.acomp(10,4),bandwidth=4)$p.value)
hist(dd)


dd <- replicate(1000,acompGOF.test(runif.acomp(20,4),runif.acomp(100,4)+acomp(c(1,2,3,1)))$p.value)

hist(dd)

# test uniformity

attach("gsi") # the uniformity test is only available as an internal function
x <- runif.acomp(100,4)
gsi.acompUniformityGOF.test.test(x)

dd <- replicate(1000,gsi.acompUniformityGOF.test.test(runif.acomp(10,4))$p.value)
hist(dd)
detach("gsi")

}
}
\keyword{htest}
