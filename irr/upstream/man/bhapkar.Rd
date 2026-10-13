\name{bhapkar}
\alias{bhapkar}
\title{Bhapkar coefficient of concordance between raters}
\description{
Calculates the Bhapkar coefficient of concordance for two raters.
}
\usage{
bhapkar(ratings)
}
\arguments{
  \item{ratings}{n*2 matrix or dataframe, n subjects 2 raters.}
}
\details{
Missing data are omitted in a listwise way.
The Bhapkar (1966) test is a more powerful alternative to the Stuart-Maxwell test.
Both tests are asymptotically equivalent and will produce comparable chi-squared
values when applied a large sample of rated objects.
}
\value{
 A list with class "irrlist" containing the following components: 
 \item{$method}{a character string describing the method.}
 \item{$subjects}{the number of data objects.}
 \item{$raters}{the number of raters.}
 \item{$irr.name}{the name of the coefficient (Chisq).}
 \item{$value}{the value of the coefficient.}
 \item{$stat.name}{the name and df of the test statistic.}
 \item{$statistic}{the value of the test statistic.}
 \item{$p.value}{the probability of the test statistic.}
}
\references{
Bhapkar, V.P. (1966). A note on the equivalence of two test criteria for hypotheses in 
categorical data. Journal of the American Statistical Association, 61, 228-235.
}
\author{Matthias Gamer}
\seealso{
  \code{\link[stats]{mcnemar.test}}, 
  \code{\link[irr]{stuart.maxwell.mh}}, 
  \code{\link[irr]{rater.bias}} 
}
\examples{
data(vision)
bhapkar(vision) # Original example used from Bhapkar (1966)
}
\keyword{misc}
