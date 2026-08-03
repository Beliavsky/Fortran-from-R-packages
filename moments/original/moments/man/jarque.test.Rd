\name{jarque.test}
\alias{jarque.test}
\title{ Jarque-Bera test for normality }
\description{
  This function performs the Jarque-Bera test on the given data sample to determine
  if the data are sample drawn from a normal population.
}
\usage{
jarque.test(x)
}
\arguments{
  \item{x}{ a numeric vector of data }
}
\details{
  Under the hypothesis of normality, data should be symmetrical (i.e. skewness should
  be equal to zero) and have skewness chose to three.  The Jarque-Bera statistic is
  chi-square distributed with two degrees of freedom.
}
\value{
A list with class \code{htest} containing the following components:

    \item{statistic}{the list containing the Jarque-Bera statistic}
    \item{p.value}{the p-value for the test.}
    \item{alternative}{a character string describing the alternative hypothesis.}
    \item{method}{a character string indicating what type of test was performed.}
    \item{data.name}{name of the data argument.}

}
\references{

 Jarque, C. M., Bera, A. K. (1980) Efficient test for normality, homoscedasticity and serial independence
 of residuals, Economic Letters, Vol. 6 Issue 3, 255-259.
 
 }
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\examples{
set.seed( 1234 )
x <- rnorm( 1000 )
jarque.test( x )
}
\keyword{ math }
\keyword{ htest }
