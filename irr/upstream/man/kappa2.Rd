\name{kappa2}
\alias{kappa2}
\title{Cohen's Kappa and weighted Kappa for two raters}
\description{
Calculates Cohen's Kappa and weighted Kappa as an index of interrater agreement between 2 raters on categorical (or ordinal) data. Own weights for the various degrees of disagreement could be specified.
}
\usage{
kappa2(ratings, weight = c("unweighted", "equal", "squared"), sort.levels = FALSE)
}
\arguments{
  \item{ratings}{n*2 matrix or dataframe, n subjects 2 raters.}
  \item{weight}{either a character string specifying one predifined set of weights or a numeric vector with own weights (see details).}
  \item{sort.levels}{boolean value describing whether factor levels should be (re-)sorted during the calculation.}
}
\details{
Missing data are omitted in a listwise way.\cr
During computation, ratings are converted to factors. Therefore, the categories are ordered accordingly. When ratings are numeric, a sorting of factor levels occurs automatically. Otherwise, levels are sorted when the function is called with sort.levels=TRUE. \cr
\code{kappa2} allows for calculating weighted Kappa coefficients. Beneath '"unweighted"' (default), predifined sets of weights are '"equal"' (all levels disagreement between raters are weighted equally) and '"squared"' (disagreements are weighted according to their squared distance from perfect agreement). The weighted Kappa coefficient with '"squared"' weights equals the product moment correlation under certain conditions.
Own weights could be specified by supplying the function with a numeric vector of weights, starting from perfect agreement to worst disagreement. The length of this vector must equal the number of rating categories.
}
\value{
A list with class '"irrlist"' containing the following components:
  \item{$method}{a character string describing the method and the weights applied for the computation of weighted Kappa.}
  \item{$subjects}{the number of subjects examined.}
  \item{$raters}{the number of raters (=2).}
  \item{$irr.name}{a character string specifying the name of the coefficient.}
  \item{$value}{value of Kappa.}
  \item{$stat.name}{a character string specifying the name of the corresponding test statistic.}
  \item{$statistic}{the value of the test statistic.}
  \item{$p.value}{the p-value for the test.}
}
\references{
Cohen, J. (1960). A coefficient of agreement for nominal scales. Educational and Psychological Measurement, 20, 37-46.\cr\cr
Cohen, J. (1968). Weighted kappa: Nominal scale agreement with provision for scaled disagreement or partial credit. Psychological Bulletin, 70, 213-220.\cr\cr
Fleiss, J.L., Cohen, J., & Everitt, B.S. (1969). Large sample standard errors of kappa and weighted kappa. Psychological Bulletin, 72, 323-327.
}
\author{Matthias Gamer}
\seealso{
  \code{\link[stats]{cor}},
  \code{\link[irr]{kappa2}},
  \code{\link[irr]{kappam.light}}
}
\examples{
data(anxiety)
kappa2(anxiety[,1:2], "squared") # predefined set of squared weights
kappa2(anxiety[,1:2], (0:5)^2)   # same result with own set of squared weights

# own weights increasing gradually with larger distance from perfect agreement
kappa2(anxiety[,1:2], c(0,1,2,4,7,11))

data(diagnoses)
# Unweighted Kappa for categorical data without a logical order
kappa2(diagnoses[,2:3])
}
\keyword{univar}
