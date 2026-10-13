\name{stuart.maxwell.mh}
\alias{stuart.maxwell.mh}
\title{Stuart-Maxwell coefficient of concordance for two raters}
\description{
 Calculates the Stuart-Maxwell coefficient of concordance for two raters.
}
\usage{
stuart.maxwell.mh(x)
}
\arguments{
 \item{x}{c x c classification matrix or matrix of classification scores
  into c categories.}
}
\details{
 \code{stuart.maxwell.mh} calculates a reliability coefficient for two raters
  classifying n objects into any number of categories. It will accept either 
  a c x c classification matrix of counts of objects falling into c categories
  or a c x n or n x c matrix of classification scores.
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
 Stuart, A.A. (1955). A test for homogeneity of the marginal distributions in a two-way classification. Biometrika, 42, 412-416.\cr\cr
 Maxwell, A.E. (1970) Comparing the classification of subjects by two independent judges. British Journal of Psychiatry, 116, 651-655.
}
\author{Jim Lemon}
\seealso{
  \code{\link[irr]{bhapkar}}, 
  \code{\link[irr]{rater.bias}}
}
\examples{
 # fake a 2xn matrix of three way classification scores
 ratings<-matrix(sample(1:3,60,TRUE), nrow=2)
 stuart.maxwell.mh(ratings)
 
 # Example used from Stuart (1955)
 data(vision)
 stuart.maxwell.mh(vision)
}
\keyword{misc}
