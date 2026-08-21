\name{rlnorm}
\alias{rlnorm.rplus}
\alias{dlnorm.rplus}
\title{The multivariate lognormal distribution}
\description{
  Generates random amounts with a multivariate lognormal distribution, or gives the 
  density of that distribution at a given point. 
}
\usage{
rlnorm.rplus(n,meanlog,varlog)
dlnorm.rplus(x,meanlog,varlog)
          }
\arguments{
  \item{n}{number of datasets to be simulated}
  \item{meanlog}{the mean-vector of the logs}
  \item{varlog}{the variance/covariance matrix of the logs}
  \item{x}{vectors in the sample space}
}
\value{
  \code{rlnorm.rplus} gives a generated random dataset of class
\code{"rplus"} following a
  lognormal distribution with logs having mean \code{meanlog} and
  variance \code{varlog}.
  \cr
  \code{dlnorm.rplus} gives the density of the distribution with respect
  to the Lesbesgue measure on R+ as a subset of R. 
}
\seealso{
  \code{\link{rnorm.acomp}}
}
\note{The main difference between \code{rlnorm.rplus} and
  \code{rnorm.aplus}
  is that rlnorm.rplus needs a logged mean. The additional difference
  for the calculation of the density by \code{dlnorm.rplus} and
  \code{dnorm.aplus} is the reference measure (a log-Lebesgue one in the
  second case).
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}, Raimon
Tolosana-Delgado}
\references{
  Aitchison, J. (1986) \emph{The Statistical Analysis of Compositional
    Data} Monographs on Statistics and Applied Probability. Chapman &
  Hall Ltd., London (UK). 416p.\cr
}
\examples{
MyVar <- matrix(c(
0.2,0.1,0.0,
0.1,0.2,0.0,
0.0,0.0,0.2),byrow=TRUE,nrow=3)
MyMean <- c(1,1,2)

plot(rlnorm.rplus(100,log(MyMean),MyVar))
plot(rnorm.aplus(100,MyMean,MyVar))
x <- rnorm.aplus(5,MyMean,MyVar)
dnorm.aplus(x,MyMean,MyVar)
dlnorm.rplus(x,log(MyMean),MyVar)

}
\keyword{multivariate}
