\name{LdaClassic}
\alias{LdaClassic}
\alias{LdaClassic.formula}
\alias{LdaClassic.default}
\title{ Linear Discriminant Analysis }
\description{
  Performs a linear discriminant analysis and returns the results as an object of class \code{LdaClassic} (aka constructor).
}
\usage{
LdaClassic(x, ...)

%%\method{LdaClassic}{formula}(formula, data = NULL, subset, na.action, \dots)
\method{LdaClassic}{default}(x, grouping, prior = proportions, tol = 1.0e-4, ...)
}
\arguments{
  \item{x}{a matrix or data frame containing the explanatory variables (training set). }
  \item{grouping}{grouping variable:  a factor specifying the class for each observation.}
  \item{prior}{prior probabilities, default to the class proportions for the training set.}
  \item{tol}{tolerance}
  \item{\dots}{arguments passed to or from other methods.}
}
%\details{}
\value{
  Returns an S4 object of class \code{LdaClassic}
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at}
}
\seealso{
     \code{\link{Lda-class}}, \code{\link{LdaClassic-class}},
}
\examples{
## Example anorexia
library(MASS)
data(anorexia)

## rrcov: LdaClassic()
lda <- LdaClassic(Treat~., data=anorexia)
predict(lda)@classification

## MASS: lda()
lda.MASS <- lda(Treat~., data=anorexia)
predict(lda.MASS)$class

## Compare the prediction results of MASS:::lda() and LdaClassic()
all.equal(predict(lda)@classification, predict(lda.MASS)$class)
}
\keyword{robust}
\keyword{multivariate}
