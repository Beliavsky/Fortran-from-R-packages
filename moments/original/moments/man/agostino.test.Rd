\name{agostino.test}
\alias{agostino.test}
\title{ D'Agostino test of skewness }
\description{
  Performs D'Agostino test for skewness in normally distributed data.
}
\usage{
agostino.test(x, alternative = c("two.sided", "less", "greater"))
}
\arguments{
  \item{x}{ a numeric vector of data values. }
  \item{alternative}{ a character string specifying the alternative hypothesis,
   must be one of '"two.sided"' (default), '"greater"' or
   '"less"'.  You can specify just the initial letter.
 }
}
\details{
  Under the hypothesis of normality, data should be symmetrical (i.e. skewness should
  be equal to zero). This test has such null hypothesis and is useful to detect a
  significant skewness in normally distributed data.
}
\value{

A list with class \code{htest} containing the following components:

	\item{statistic}{the list containing skewness estimator and its transformation.}
	\item{p.value}{the p-value for the test.}
	\item{alternative}{a character string describing the alternative hypothesis.}
	\item{method}{a character string indicating what type of test was performed.}
	\item{data.name}{name of the data argument.}

}

\references{ 

D'Agostino, R.B. (1970). Transformation to Normality of the Null Distribution of G1. 
 Biometrika, 57, 3, 679-681.

}
\author{ Lukasz Komsta }

\seealso{ \code{\link{skewness}} }
\examples{
set.seed(1234)
x = rnorm(1000)
skewness(x)
agostino.test(x)
}
\keyword{htest}
