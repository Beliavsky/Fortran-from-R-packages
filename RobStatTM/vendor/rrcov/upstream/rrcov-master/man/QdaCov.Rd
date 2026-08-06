\name{QdaCov}
\alias{QdaCov}
\alias{QdaCov.formula}
\alias{QdaCov.default}
\title{ Robust Quadratic Discriminant Analysis }
\description{
    Performs robust quadratic discriminant analysis and returns 
    the results as an object of class \code{QdaCov} (aka constructor). 
}
\usage{
QdaCov(x, ...)

\method{QdaCov}{default}(x, grouping, prior = proportions, tol = 1.0e-4,
                 method = CovControlMcd(), ...)
}
\arguments{
  \item{x}{a matrix or data frame containing the explanatory variables (training set). }
  \item{grouping}{grouping variable:  a factor specifying the class for each observation.}
  \item{prior}{prior probabilities, default to the class proportions for the training set.}
  \item{tol}{tolerance}
  \item{method}{method}
  \item{\dots}{arguments passed to or from other methods}
}
\details{
  details
}
\value{
  Returns an S4 object of class \code{QdaCov}
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} 
}
\section{Warning }{Still an experimental version!}
\seealso{  \code{\link[rrcov]{CovMcd}} }
\examples{
## Example anorexia
library(MASS)
data(anorexia)

## start with the classical estimates
qda <- QdaClassic(Treat~., data=anorexia)
predict(qda)@classification

## try now the robust LDA with the default method (MCD with pooled whitin cov matrix)
rqda <- QdaCov(Treat~., data= anorexia)
predict(rqda)@classification

## try the other methods
QdaCov(Treat~., data= anorexia, method="sde")
QdaCov(Treat~., data= anorexia, method="M")
QdaCov(Treat~., data= anorexia, method=CovControlOgk())

}
\keyword{robust}
\keyword{multivariate}
