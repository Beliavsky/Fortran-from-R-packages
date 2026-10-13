\name{N2.cohen.kappa}
\alias{N2.cohen.kappa}
\title{Sample Size Calculation for Cohen's Kappa Statistic with more than one category}
\description{
 This function calculates the required sample size for the Cohen's Kappa
 statistic when two raters have the same marginal. 
 Note that any value of "kappa under null" in the interval [-1,1] is
 acceptable (i.e. k0=0 is a valid null hypothesis).
}
\usage{
 N2.cohen.kappa(mrg, k1, k0, alpha=0.05, power=0.8, twosided=FALSE)
}

\arguments{
  \item{mrg}{a vector of marginal probabilities given by raters}
  \item{k1}{the true Cohen's Kappa statistic}
  \item{k0}{the value of kappa under the null hypothesis}
  \item{alpha}{type I error of test}
  \item{power}{the desired power to detect the difference between
   true kappa and hypothetical kappa}
  \item{twosided}{TRUE if test is two-sided}
}

\value{
    Returns required sample size.
}
\references{
  Flack, V.F., Afifi, A.A., Lachenbruch, P.A., & Schouten, H.J.A. (1988). Sample size determinations for the two rater kappa statistic. Psychometrika, 53, 321-325.
}
\author{Puspendra Singh and Jim Lemon}
\seealso{
  \code{\link[irr]{N.cohen.kappa}},
  \code{\link[irr]{kappa2}}
}
\examples{
  # Testing H0: kappa = 0.4 vs. HA: kappa > 0.4 (=0.6) given that
  # Marginal Probabilities by two raters are (0.2, 0.25, 0.55).
  #
  # one sided test with 80% power:
  N2.cohen.kappa(c(0.2, 0.25, 0.55), k1=0.6, k0=0.4)
  
  # one sided test with 90% power:
  N2.cohen.kappa(c(0.2, 0.25, 0.55), k1=0.6, k0=0.4, power=0.9)	  

  # Marginal Probabilities by two raters are (0.2, 0.05, 0.2, 0.05, 0.2, 0.3)
  # Testing H0: kappa = 0.1 vs. HA: kappa > 0.1 (=0.5) given that
  #
  # one sided test with 80% power:
  N2.cohen.kappa(c(0.2, 0.05, 0.2, 0.05, 0.2, 0.3), k1=0.5, k0=0.1)				
}
\keyword{misc}
