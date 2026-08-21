\name{ccompgof}
\alias{PoissonGOF.test}
\alias{ccompPoissonGOF.test}
\alias{ccompMultinomialGOF.test}
\title{Compositional Goodness of fit test}
\description{
  Goodness of fit tests for count compositional data.
}
\usage{
PoissonGOF.test(x,lambda=mean(x),R=999,estimated=missing(lambda))
ccompPoissonGOF.test(x,simulate.p.value=TRUE,R=1999)
ccompMultinomialGOF.test(x,simulate.p.value=TRUE,R=1999)
}
\arguments{
  \item{x}{a dataset integer numbers (PoissonGOF) or count compositions (compPoissonGOF)}
  \item{lambda}{the expected value to check against}
  \item{R}{The number of replicates to compute the distribution
    of the test statistic}
  \item{estimated}{states whether the lambda parameter should be
    considered as estimated for the computation of the p-value.}
  \item{simulate.p.value}{should all p-values be infered by simulation.}
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

erg <- acompGOF.test(x,y)
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


x <- runif.acomp(100,4)
acompUniformityGOF.test(x)

dd <- replicate(1000,acompUniformityGOF.test(runif.acomp(10,4))$p.value)

hist(dd)

}
}
\keyword{htest}
