\name{gausstest}
\alias{Gauss.test}
\title{Classical Gauss Test}
\description{
  One and two sample Gauss test for equal mean of normal random variates
  with known variance.
}
\usage{
Gauss.test(x,y=NULL,mean=0,sd=1,alternative = c("two.sided", "less", "greater"))
          }
\arguments{
  \item{x}{a numeric vector providing the first dataset }
  \item{y}{optional second dataset}
  \item{mean}{the mean to compare with}
  \item{sd}{the known standard deviation }
  \item{alternative}{the alternative to be used in the test}
}
\value{
  A classical \code{"htest"} object
  \item{data.name}{The name of the dataset as specified}
  \item{method}{a name for the test used}
  \item{parameter}{the mean and variance provided to the test}
  \item{alternative}{an empty string}
  \item{p.value}{The p.value computed for this test}
}
\details{
The Gauss test is in every Text-Book, but not in R, because it is nearly
never used. However it is included here for educational purposes.  
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}

\seealso{
  \code{\link{t.test}} 
}
\examples{
x <- rnorm(100)
y <- rnorm(100)
Gauss.test(x,y)

}
\keyword{htest}
