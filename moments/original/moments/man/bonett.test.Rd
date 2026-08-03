\name{bonett.test}
\alias{bonett.test}
\title{ Bonett-Seier test of Geary's kurtosis }
\description{
 This function performs Bonett-Seier test of Geary's measure of kurtosis
 for normally distributed data.
}
\usage{
bonett.test(x, alternative = c("two.sided", "less", "greater"))
}
\arguments{
  \item{x}{ a numeric vector of data values. }
  \item{alternative}{ a character string specifying the alternative hypothesis,
   must be one of '"two.sided"' (default), '"greater"' or
   '"less"'. You can specify just the initial letter.
 }
}
\details{
  Under the hypothesis of normality, data should have Geary's kurtosis 
  equal to \code{sqrt(2/pi)} (0.7979). This test has such null hypothesis and is 
  useful to detect a significant difference of Geary's kurtosis in normally 
  distributed data.
}
\value{

A list with class \code{htest} containing the following components:

	\item{statistic}{the list containing Geary's kurtosis estimator and its transformation.}
	\item{p.value}{the p-value for the test.}
	\item{alternative}{a character string describing the alternative hypothesis.}
	\item{method}{a character string indicating what type of test was performed.}
	\item{data.name}{name of the data argument.}

}

\references{

Bonett, D.G., Seier, E. (2002) A test of normality with high uniform power. Computational
 Statistics and Data Analysis, 40, 435-445.

}
\author{ Lukasz Komsta }
\seealso{ \code{\link{geary}}}
\examples{
set.seed(1234)
x = rnorm(1000)
geary(x)
bonett.test(x)
}
\keyword{htest}

