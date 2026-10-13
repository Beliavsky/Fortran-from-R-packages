\name{meanrho}
\alias{meanrho}
\title{Mean of bivariate rank correlations between raters}
\description{
Computes the mean of bivariate Spearman's rho rank correlations between raters as an index of the interrater reliability of ordinal data.
}
\usage{
meanrho(ratings, fisher = TRUE)
}
\arguments{
  \item{ratings}{n*m matrix or dataframe, n subjects m raters.}
  \item{fisher}{a logical indicating whether the correlation coefficients should be Fisher z-standardized before averaging.}
}
\details{
Missing data are omitted in a listwise way.\cr
The mean of bivariate rank correlations should not be used as an index of interrater reliability when ties within raters occur.\cr
The null hypothesis r=0 could only be tested when Fisher z-standardized values are used for the averaging.\cr
When computing Fisher z-standardized values, perfect correlations are omitted before averaging because z equals +/-Inf in that case.
}
\value{
A list with class '"irrlist"' containing the following components:
  \item{$method}{a character string describing the method applied for the computation of interrater reliability.}
  \item{$subjects}{the number of subjects examined.}
  \item{$raters}{the number of raters.}
  \item{$irr.name}{a character string specifying the name of the coefficient.}
  \item{$value}{coefficient of interrater reliability.}
  \item{$stat.name}{a character string specifying the name of the corresponding test statistic.}
  \item{$statistic}{the value of the test statistic.}
  \item{$p.value}{the p-value for the test.}
  \item{$error}{a character specifying whether correlations were dropped before the computation of the Fisher z-standardized average. Additionally, a warning message is created if ties were found within raters.}
}
\author{Matthias Gamer}
\seealso{
  \code{\link[stats]{cor}},
  \code{\link[irr]{kendall}}
}
\examples{
data(anxiety)
meanrho(anxiety, TRUE)
}
\keyword{univar}
