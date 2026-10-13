\name{cochran.test}
\alias{cochran.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Test for outlying or inlying variance }
\description{

This test is useful to check if largest variance in several groups of
data is "outlying" and this group should be rejected. Alternatively,
if one group has very small variance, we can test for "inlying" variance.

}
\usage{
cochran.test(object, data, inlying = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ A vector of variances or formula. }
  \item{data}{ If object is a vector, data should be another vector, giving number of data in
each corresponding group. If object is a formula, data should be a dataframe.}
  \item{inlying}{ Test smallest variance instead of largest.}
}
\details{

The corresponding p-value is calculated using \code{pcochran} function.

}
\value{

A list with class \code{htest} containing the following components:

	\item{statistic}{the value of Cochran-statistic.}
	\item{p.value}{the p-value for the test.}
	\item{alternative}{a character string describing the alternative hypothesis.}
	\item{method}{a character string indicating what type of test was performed.}
	\item{data.name}{name of the data argument.}
	\item{estimate}{vector of variance estimates}


}
\references{
Snedecor, G.W., Cochran, W.G. (1980). Statistical Methods (seventh edition).
		 Iowa State University Press, Ames, Iowa.
 }
\author{ Lukasz Komsta }

\seealso{ \code{\link{qcochran}} }
\examples{
set.seed(1234)
x=rnorm(100)
d=data.frame(x=x,group=rep(1:10,10))
cochran.test(x~group,d)
cochran.test(x~group,d,inlying=TRUE)
x=runif(5)
cochran.test(x,rep(5,5))
cochran.test(x,rep(100,5))
}
\keyword{htest}