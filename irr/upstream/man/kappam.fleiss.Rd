\name{kappam.fleiss}
\alias{kappam.fleiss}
\title{Fleiss' Kappa for m raters}
\description{
Computes Fleiss' Kappa as an index of interrater agreement between m raters on categorical data. Additionally, category-wise Kappas could be computed.
}
\usage{
kappam.fleiss(ratings, exact = FALSE, detail = FALSE)
}
\arguments{
  \item{ratings}{n*m matrix or dataframe, n subjects m raters.}
  \item{exact}{a logical indicating whether the exact Kappa (Conger, 1980) or the Kappa described by Fleiss (1971) should be computed.}
  \item{detail}{a logical indicating whether category-wise Kappas should be computed}
}
\details{
Missing data are omitted in a listwise way.\cr
The coefficient described by Fleiss (1971) does not reduce to Cohen's Kappa (unweighted) for m=2 raters. Therefore, the exact Kappa coefficient, which is slightly higher in most cases, was proposed by Conger (1980).\cr
The null hypothesis Kappa=0 could only be tested using Fleiss' formulation of Kappa.
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
  \item{$detail}{a table with category-wise kappas and the corresponding test statistics.}
}
\references{
Conger, A.J. (1980). Integration and generalisation of Kappas for multiple raters. Psychological Bulletin, 88, 322-328.\cr\cr
Fleiss, J.L. (1971). Measuring nominal scale agreement among many raters. Psychological Bulletin, 76, 378-382.\cr\cr
Fleiss, J.L., Levin, B., & Paik, M.C. (2003). Statistical Methods for Rates and Proportions, 3rd Edition. New York: John Wiley & Sons.
}
\author{Matthias Gamer}
\seealso{
  \code{\link[irr]{kappa2}},
  \code{\link[irr]{kappam.light}}
}
\examples{
data(diagnoses)
kappam.fleiss(diagnoses)               # Fleiss' Kappa
kappam.fleiss(diagnoses, exact=TRUE)   # Exact Kappa
kappam.fleiss(diagnoses, detail=TRUE)  # Fleiss' and category-wise Kappa

kappam.fleiss(diagnoses[,1:4])         # Fleiss' Kappa of raters 1 to 4
}
\keyword{univar}
