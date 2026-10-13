\name{finn}
\alias{finn}
\title{Finn coefficient for oneway and twoway models}
\description{
Computes the Finn coefficient as an index of the interrater reliability of quantitative data. Additionally, F-test and confidence interval are computed.
}
\usage{
finn(ratings, s.levels, model = c("oneway", "twoway"))
}
\arguments{
  \item{ratings}{n*m matrix or dataframe, n subjects m raters.}
  \item{s.levels}{the number of different rating categories.}
  \item{model}{a character string specifying if a '"oneway"' model (default) with row effects random, or a '"twoway"' model with column and row effects random should be applied. You can specify just the initial letter.}
}
\details{
Missing data are omitted in a listwise way.\cr
The Finn coefficient is especially useful, when variance between raters is low (i.e. agreement is high).\cr
For the computation it could be specified if only the subjects are considered as random effects ('"oneway"' model) or if subjects and raters are randomly chosen from a bigger pool of persons ('"twoway"' model).
}
\value{
A list with class '"irrlist"' containing the following components:
  \item{$method}{a character string describing the method applied for the computation of interrater reliability.}
  \item{$subjects}{the number of subjects examined.}
  \item{$raters}{the number of raters.}
  \item{$irr.name}{a character string specifying the name of the coefficient.}
  \item{$value}{coefficient of interrater reliability.}
  \item{$stat.name}{a character string specifying the name and the df of the corresponding F-statistic.}
  \item{$statistic}{the value of the test statistic.}
  \item{$p.value}{the p-value for the test.}
}
\references{
Finn, R.H. (1970). A note on estimating the reliability of categorical data. Educational and Psychological Measurement, 30, 71-76.
}
\author{Matthias Gamer}
\seealso{
  \code{\link[irr]{icc}}, 
  \code{\link[irr]{meancor}}, 
  \code{\link[irr]{robinson}}
}
\examples{
data(video)
finn(video, 6, model="twoway")
}
\keyword{univar}
