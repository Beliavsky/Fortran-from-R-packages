\name{grubbs.test}
\alias{grubbs.test}
\title{ Grubbs tests for one or two outliers in data sample}
\description{

Performs Grubbs' test for one outlier, two outliers on one tail,
or two outliers on opposite tails, in small sample.

}
\usage{
grubbs.test(x, type = 10, opposite = FALSE, two.sided = FALSE)
}
\arguments{
  \item{x}{ a numeric vector for data values. }
  \item{opposite}{  a logical indicating whether you want to check not the value with 
largest difference from the mean, but opposite (lowest, if most suspicious is highest etc.)}
  \item{type}{ Integer value indicating test variant. 10 is a test for one outlier (side is
	detected automatically and can be reversed by \code{opposite} parameter). 11 is a test
	for two outliers on opposite tails, 20 is test for two outliers in one tail.  }
  \item{two.sided}{ Logical value indicating if there is a need to treat this test as two-sided. }
}
\details{

The function can perform three tests given and discussed by Grubbs (1950).

First test (10) is used to detect if the sample dataset contains one outlier, statistically different than the other values.
Test is based by calculating score of this outlier G (outlier minus mean and divided by sd) and
comparing it to appropriate critical values. Alternative method is calculating ratio of variances
of two datasets - full dataset and dataset without outlier. The obtained value called U is bound
with G by simple formula.

Second test (11) is used to check if lowest and highest value are two outliers on opposite tails
of sample. It is based on calculation of ratio of range to standard deviation of the sample. 

Third test (20) calculates ratio of variance of full sample and sample without two extreme observations.
It is used to detect if dataset contains two outliers on the same tail.

The p-values are calculated using \code{\link{qgrubbs}} function.

}
\value{

	\item{statistic}{the value statistic. For type 10 it is difference between outlier and
	the mean divided by standard deviation, and for type 20 it is sample range divided
	by standard deviation. Additional value U is ratio of sample variances with and withour suspicious outlier.
	According to Grubbs (1950) these values for type 10 are bound by simple formula and only one of them
	can be used, but function gives both. For type 20 the G is the same as U.}
	\item{p.value}{the p-value for the test.}
	\item{alternative}{a character string describing the alternative hypothesis.}
	\item{method}{a character string indicating what type of test was performed.}
	\item{data.name}{name of the data argument.}

}
\references{ 

Grubbs, F.E. (1950). Sample Criteria for testing outlying observations. Ann. Math. Stat. 21, 1, 27-58.

}
\author{ Lukasz Komsta }

\seealso{ \code{\link{dixon.test}}, \code{\link{chisq.out.test}} }

\examples{

set.seed(1234)
x = rnorm(10)
grubbs.test(x)
grubbs.test(x,type=20)
grubbs.test(x,type=11)

}
\keyword{htest}

