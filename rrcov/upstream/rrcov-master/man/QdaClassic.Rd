\name{QdaClassic}
\Rdversion{1.1}
\alias{QdaClassic}
\alias{QdaClassic.formula}
\alias{QdaClassic.default}
\title{ Quadratic Discriminant Analysis }
\description{
  Performs quadratic discriminant analysis and returns the results as an object 
  of class \code{QdaClassic} (aka constructor).
}
\usage{
QdaClassic(x, ...)

%%\method{QdaClassic}{formula}(formula, data = NULL, subset, na.action, \dots)
\method{QdaClassic}{default}(x, grouping, prior = proportions, tol = 1.0e-4, ...)
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
  Returns an S4 object of class \code{QdaClassic}
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} 
}
\seealso{
     \code{\link{Qda-class}}, \code{\link{QdaClassic-class}}, 
}
%\examples{}
\keyword{robust}
\keyword{multivariate}
