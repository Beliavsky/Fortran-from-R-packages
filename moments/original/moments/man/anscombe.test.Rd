\name{anscombe.test}
\alias{anscombe.test}
\title{ Anscombe-Glynn test of kurtosis }
\description{
  Performs Anscombe-Glynn test of kurtosis for normal samples
}
\usage{
anscombe.test(x, alternative = c("two.sided", "less", "greater"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a numeric vector of data values. }
  \item{alternative}{ a character string specifying the alternative hypothesis,
   must be one of '"two.sided"' (default), '"greater"' or
   '"less"'. You can specify just the initial letter.
 }
}
\details{
  Under the hypothesis of normality, data should have kurtosis equal to 3. 
  This test has such null hypothesis and is useful to detect a
  significant difference of kurtosis in normally distributed data.
}
\value{

A list with class \code{htest} containing the following components:

	\item{statistic}{the list containing kurtosis estimator and its transformation.}
	\item{p.value}{the p-value for the test.}
	\item{alternative}{a character string describing the alternative hypothesis.}
	\item{method}{a character string indicating what type of test was performed.}
	\item{data.name}{name of the data argument.}

}

\references{

Anscombe, F.J., Glynn, W.J. (1983) Distribution of kurtosis statistic for normal statistics.
 Biometrika, 70, 1, 227-234

}
\author{ Lukasz Komsta }
\seealso{ \code{\link{kurtosis}}}
\examples{
set.seed(1234)
x = rnorm(1000)
kurtosis(x)
anscombe.test(x)
}
\keyword{htest}

