\name{PcaGrid}
\alias{PcaGrid}
\alias{PcaGrid.formula}
\alias{PcaGrid.default}
\title{ Robust Principal Components based on Projection Pursuit (PP): GRID search Algorithm }
\description{
    Computes an approximation of the PP-estimators for PCA using the grid search algorithm in the plane.
}
\usage{
    PcaGrid(x, ...)
    \method{PcaGrid}{default}(x, k = 0, kmax = ncol(x), 
        scale=FALSE, na.action = na.fail, crit.pca.distances = 0.975, trace=FALSE, \dots)
    \method{PcaGrid}{formula}(formula, data = NULL, subset, na.action, \dots)
}
\arguments{
  \item{formula}{a formula with no response variable, referring only to
    numeric variables.}
  \item{data}{an optional data frame (or similar: see
    \code{\link{model.frame}}) containing the variables in the
    formula \code{formula}.}
  \item{subset}{an optional vector used to select rows (observations) of the
    data matrix \code{x}.}
  \item{na.action}{a function which indicates what should happen
    when the data contain \code{NA}s.  The default is set by
    the \code{na.action} setting of \code{\link{options}}, and is
    \code{\link{na.fail}} if that is unset. The default is \code{\link{na.omit}}.}
  \item{\dots}{arguments passed to or from other methods.}
  \item{x}{a numeric matrix (or data frame) which provides
    the data for the principal components analysis.}
  \item{k}{number of principal components to compute. If \code{k} is missing, 
    or \code{k = 0}, it is set to the number of columns of the data. 
    It is preferable to investigate the scree plot in order to choose the number 
    of components and then run again. Default is \code{k=0}. }
  \item{kmax}{maximal number of principal components to compute.
    Default is \code{kmax=10}. If \code{k} is provided, \code{kmax} 
    does not need to be specified, unless \code{k} is larger than 10.}
  \item{scale}{a value indicating whether and how the variables should be 
  scaled. If \code{scale = FALSE} (default) or \code{scale = NULL} no scaling is 
  performed (a vector of 1s is returned in the \code{scale} slot). 
  If \code{scale = TRUE} the data are scaled to have unit variance. Alternatively it can 
  be a function like \code{sd} or \code{mad} or a vector of length equal 
  the number of columns of \code{x}. The value is passed to the underlying function
  and the result returned is stored in the \code{scale} slot. 
  Default is \code{scale = FALSE}}
  \item{crit.pca.distances}{criterion to use for computing the cutoff values 
    for the orthogonal and score distances. Default is 0.975.}
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}}
}
\details{
\code{PcaGrid}, serving as a constructor for objects of class \code{\link{PcaGrid-class}} 
is a generic function with "formula" and "default" methods. For details see \code{\link[pcaPP]{PCAgrid}} and the relevant references.
}
\value{
  An S4 object of class \code{\link{PcaGrid-class}} which is a subclass of the 
  virtual class \code{\link{PcaRobust-class}}. 
}
\references{
  C. Croux, P. Filzmoser, M. Oliveira, (2007).
  Algorithms for Projection-Pursuit Robust Principal Component Analysis,
  \emph{Chemometrics and Intelligent Laboratory Systems}, 87, 225.

  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} 
}
\examples{
    # multivariate data with outliers
    library(mvtnorm)
    x <- rbind(rmvnorm(200, rep(0, 6), diag(c(5, rep(1,5)))),
                rmvnorm( 15, c(0, rep(20, 5)), diag(rep(1, 6))))
    # Here we calculate the principal components with PCAgrid
    pc <- PcaGrid(x, 6)
    # we could draw a biplot too:
    biplot(pc)
    
    # we could use another objective function, and 
    # maybe only calculate the first three principal components:
    pc <- PcaGrid(x, 3, method="qn")
    biplot(pc)
    
    # now we want to compare the results with the non-robust principal components
    pc <- PcaClassic(x, k=3)
    # again, a biplot for comparision:
    biplot(pc)
}
\keyword{robust}
\keyword{multivariate}
