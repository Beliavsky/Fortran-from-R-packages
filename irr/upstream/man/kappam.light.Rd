\name{kappam.light}
\alias{kappam.light}
\title{Light's Kappa for m raters}
\description{
Computes Light's Kappa as an index of interrater agreement between m raters on categorical data.
}
\usage{
kappam.light(ratings)
}
\arguments{
  \item{ratings}{n*m matrix or dataframe, n subjects m raters.}
}
\details{
Missing data are omitted in a listwise way.\cr
Light's Kappa equals the average of all possible combinations of bivariate Kappas between raters.
}
\value{
A list with class '"irrlist"' containing the following components:
  \item{$method}{a character string describing the method applied for the computation of interrater reliability.}
  \item{$subjects}{the number of subjects examined.}
  \item{$raters}{the number of raters.}
  \item{$irr.name}{a character string specifying the name of the coefficient.}
  \item{$value}{value of Kappa.}
  \item{$stat.name}{a character string specifying the name of the corresponding test statistic.}
  \item{$statistic}{the value of the test statistic.}
  \item{$p.value}{the p-value for the test.}
}
\references{
Conger, A.J. (1980). Integration and generalisation of Kappas for multiple raters. Psychological Bulletin, 88, 322-328.\cr\cr
Light, R.J. (1971). Measures of response agreement for qualitative data: Some generalizations and alternatives. Psychological Bulletin, 76, 365-377.
}
\author{Matthias Gamer}
\seealso{
  \code{\link[irr]{kappa2}},
  \code{\link[irr]{kappam.fleiss}}
}
\examples{
data(diagnoses)
kappam.light(diagnoses)   # Light's Kappa
}
\keyword{univar}
