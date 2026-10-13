\name{kendall}
\alias{kendall}
\title{Kendall's coefficient of concordance W}
\description{
Computes Kendall's coefficient of concordance as an index of interrater reliability of ordinal data. The coefficient could be corrected for ties within raters.
}
\usage{
kendall(ratings, correct = FALSE)
}
\arguments{
  \item{ratings}{n*m matrix or dataframe, n subjects m raters.}
  \item{correct}{a logical indicating whether the coefficient should be corrected for ties within raters.}
}
\details{
Missing data are omitted in a listwise way.\cr
Kendall's W should be corrected for ties if raters did not use a true ranking order for the subjects.\cr
A test for the significance of Kendall's W is only valid for large samples. 
}
\value{
A list with class '"irrlist"' containing the following components:
  \item{$method}{a character string describing the method applied for the computation of interrater reliability.}
  \item{$subjects}{the number of subjects examined.}
  \item{$raters}{the number of raters.}
  \item{$irr.name}{a character string specifying the name of the coefficient.}
  \item{$value}{coefficient of interrater reliability.}
  \item{$stat.name}{a character string specifying the name and the df of the corresponding chi-squared test.}
  \item{$statistic}{the value of the test statistic.}
  \item{$p.value}{the p-value for the test.}
  \item{$error}{the character string of a warning message if ties were found within raters.}
}
\references{
Kendall, M.G. (1948). Rank correlation methods. London: Griffin.
}
\author{Matthias Gamer}
\seealso{
  \code{\link[stats]{cor}},
  \code{\link[irr]{meanrho}}
}
\examples{
data(anxiety)
kendall(anxiety, TRUE)
}
\keyword{univar}
