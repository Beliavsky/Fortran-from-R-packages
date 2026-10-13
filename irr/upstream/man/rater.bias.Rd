\name{rater.bias}
\alias{rater.bias}
\title{Coefficient of rater bias}
\description{
 Calculates a coefficient of systematic bias between two raters.
}
\usage{
rater.bias(x)
}
\arguments{
 \item{x}{c x c classification matrix or 2 x n or n x 2 matrix of classification
  scores into c categories.}
}
\details{
 \code{rater.bias} calculates a reliability coefficient for two raters
  classifying n objects into any number of categories. It will accept either 
  a c x c classification matrix of counts of objects falling into c categories
  or a 2 x n or n x 2 matrix of classification scores.\cr
  The function returns the absolute value of the triangular off-diagnonal 
  sum ratio of the cxc classification table and the corresponding test statistic.
  A systematic bias between raters can be assumed when the ratio substantially 
  deviates from 0.5 while yielding a significant Chi-squared statistic.
}
\value{
 \item{method}{Name of the method}
 \item{subjects}{Number of subjects}
 \item{raters}{Number of raters (2)}
 \item{irr.name}{Name of the coefficient: ratio of triangular off-diagnonal sums}
 \item{value}{Value of the coefficient}
 \item{stat.name}{Name of the test statistic}
 \item{statistic}{Value of the test statistic}
 \item{p.value}{the probability of the df 1 Chi-square variable}
}
\references{
 Bishop Y.M.M., Fienberg S.E., & Holland P.W. (1978). Discrete multivariate analysis: theory and practice. Cambridge, Massachusetts: MIT Press.
}
\author{Jim Lemon}
\seealso{
  \code{\link[stats]{mcnemar.test}}
}
\examples{
 # fake a 2xn matrix of three way classification scores
 ratings <- matrix(sample(1:3,60,TRUE), nrow=2)
 rater.bias(ratings)

 # Example from Bishop, Fienberg & Holland (1978), Table 8.2-1
 data(vision)
 rater.bias(vision)                      
}
\keyword{misc}
