\name{chisq.out.test}
\alias{chisq.out.test}
\title{ Chi-squared test for outlier }
\description{
  Performs a chisquared test for detection of one outlier in a vector.
}
\usage{
chisq.out.test(x, variance=var(x), opposite = FALSE)
}
\arguments{
  \item{x}{ a numeric vector for data values. }
  \item{variance}{ known variance of population. if not given, estimator from sample is taken,
but there is not so much sense in such test (it is similar to z-scores) }
  \item{opposite}{  a logical indicating whether you want to check not the value with 
largest difference from the mean, but opposite (lowest, if most suspicious is highest etc.)}
}
\details{

This function performs a simple test for one outlier, based on chisquared distribution of squared
differences between data and sample mean. It assumes known variance of population. 
It is rather not recommended today for routine use, because several more powerful 
tests are implemented (see other functions mentioned below). 
It was discussed by Dixon (1950) for the first time, as one of the tests taken into account by him.

}
\value{

A list with class \code{htest} containing the following components:

	\item{statistic}{the value of chisquared-statistic.}
	\item{p.value}{the p-value for the test.}
	\item{alternative}{a character string describing the alternative hypothesis.}
	\item{method}{a character string indicating what type of test was performed.}
	\item{data.name}{name of the data argument.}

}
\references{ 

Dixon, W.J. (1950). Analysis of extreme values. Ann. Math. Stat. 21, 4, 488-506.

}
\author{ Lukasz Komsta }
\note{ 

This test is known to reject only extreme outliers, if no known variance is specified.

}

\seealso{ \code{\link{dixon.test}}, \code{\link{grubbs.test}} }
\examples{

set.seed(1234)
x = rnorm(10)
chisq.out.test(x)
chisq.out.test(x,opposite=TRUE)

}
\keyword{htest}
