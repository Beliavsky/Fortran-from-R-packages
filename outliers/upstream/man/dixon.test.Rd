\name{dixon.test}
\alias{dixon.test}
\title{ Dixon tests for outlier }
\description{
 Performs several variants of Dixon test for detecting outlier in data sample.
}
\usage{
dixon.test(x, type = 0, opposite = FALSE, two.sided = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a numeric vector for data values. }
  \item{opposite}{  a logical indicating whether you want to check not the value with 
largest difference from the mean, but opposite (lowest, if most suspicious is highest etc.)}
  \item{type}{ an integer specyfying the variant of test to be performed. Possible values are
compliant with these given by Dixon (1950): 10, 11, 12, 20, 21. If this value is set to zero,
a variant of the test is chosen according to sample size (10 for 3-7, 11 for 8-10, 21 for 11-13, 
22 for 14 and more). The lowest or highest value is selected automatically, and can be reversed
used \code{opposite} parameter.
}
  \item{two.sided}{ treat test as two-sided (default). }
}
\details{
The p-value is calculating by interpolation using \code{\link{qdixon}} and \code{\link{qtable}}.
According to Dixon (1951) conclusions, the critical values can be obtained numerically only for n=3.
Other critical values are obtained by simulations, taken from original Dixon's paper, and
regarding corrections given by Rorabacher (1991).
}
\value{

A list with class \code{htest} containing the following components:

	\item{statistic}{the value of Dixon Q-statistic.}
	\item{p.value}{the p-value for the test.}
	\item{alternative}{a character string describing the alternative hypothesis.}
	\item{method}{a character string indicating what type of test was performed.}
	\item{data.name}{name of the data argument.}

}
\references{ 

Dixon, W.J. (1950). Analysis of extreme values. Ann. Math. Stat. 21, 4, 488-506.

Dixon, W.J. (1951). Ratios involving extreme values. Ann. Math. Stat. 22, 1, 68-78.

Rorabacher, D.B. (1991). Statistical Treatment for Rejection of Deviant Values: Critical Values of Dixon Q Parameter and Related Subrange Ratios at the 95 percent Confidence Level. Anal. Chem. 83, 2, 139-146.

}
\author{ Lukasz Komsta }

\seealso{ \code{\link{chisq.out.test}}, \code{\link{grubbs.test}} }
\examples{

set.seed(1234)
x = rnorm(10)
dixon.test(x)
dixon.test(x,opposite=TRUE)
dixon.test(x,type=10)

}
\keyword{htest}
