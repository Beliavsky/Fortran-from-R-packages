\name{PcaCov}
\alias{PcaCov}
\alias{PcaCov.formula}
\alias{PcaCov.default}
\title{ Robust PCA based on a robust covariance matrix }
\description{
Robust PCA are obtained by replacing the classical covariance matrix 
by a robust covariance estimator. This can be one of the available 
in \code{rrcov} estimators, i.e. MCD, OGK, M or S estimator.
}
\usage{
PcaCov(x, ...)
\method{PcaCov}{default}(x, k = ncol(x), kmax = ncol(x), cov.control=CovControlMcd(), 
    scale = FALSE, signflip = TRUE, crit.pca.distances = 0.975, trace=FALSE, \dots)
\method{PcaCov}{formula}(formula, data = NULL, subset, na.action, \dots)
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
    or \code{k = 0}, the algorithm itself will determine the number of 
    components by finding such \code{k} that \eqn{l_k/l_1 >= 10.E-3} and 
    \eqn{\Sigma_{j=1}^k l_j/\Sigma_{j=1}^r l_j >= 0.8}. 
    It is preferable to investigate the scree plot in order to choose the number 
    of components and then run again. Default is \code{k=ncol(x)}. }
  \item{kmax}{maximal number of principal components to compute.
    Default is \code{kmax=10}. If \code{k} is provided, \code{kmax} 
    does not need to be specified, unless \code{k} is larger than 10.}
  \item{cov.control}{specifies which covariance estimator to use by providing 
  a \code{\link{CovControl-class}} object. 
  The default is \code{\link{CovControlMcd-class}} which will indirectly call \code{\link{CovMcd}}. 
  If \code{cov.control=NULL} is specified, the classical estimates will be used by calling \code{\link{CovClassic}}}.
    \item{scale}{a value indicating whether and how the variables should be scaled
    to have unit variance (only possible if there are no constant 
    variables). If \code{scale=FALSE} (default) or \code{scale=NULL} no scaling is 
    performed (a vector of 1s is returned in the scale slot). If \code{scale=TRUE}
    the data are scaled by the estimator used to compute the covariance matrix 
    (\code{MCD} by default). Alternatively it can be a function 
    like \code{sd} or \code{Qn} or a vector of length equal the number of columns 
    of \code{x}. The value is passed to the underlying function and the result 
    returned is stored in the scale slot. Default is \code{scale=FALSE}.}
    \item{signflip}{a logical value indicating wheather to try to solve the sign indeterminancy of the loadings -   
  ad hoc approach setting the maximum element in a singular vector to be positive. Default is \code{signflip = TRUE}}
  
  \item{crit.pca.distances}{criterion to use for computing the cutoff values for the orthogonal and score distances. Default is 0.975.}
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}}
}
\details{
\code{PcaCov}, serving as a constructor for objects of class \code{\link{PcaCov-class}} 
is a generic function with "formula" and "default" methods. For details see the relevant references.
}
\value{
  An S4 object of class \code{\link{PcaCov-class}} which is a subclass of the 
  virtual class \code{\link{PcaRobust-class}}. 
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} 
}
%\note{}
\examples{
## PCA of the Hawkins Bradu Kass's Artificial Data
##  using all 4 variables
    data(hbk)
    pca <- PcaCov(hbk)
    pca

## Compare with the classical PCA
    prcomp(hbk)

## or  
    PcaClassic(hbk)
    
## If you want to print the scores too, use
    print(pca, print.x=TRUE)

## Using the formula interface
    PcaCov(~., data=hbk)

## To plot the results:

    plot(pca)                    # distance plot
    pca2 <- PcaCov(hbk, k=2)  
    plot(pca2)                   # PCA diagnostic plot (or outlier map)
    
## Use the standard plots available for for prcomp and princomp
    screeplot(pca)    
    biplot(pca)    
}
\keyword{robust}
\keyword{multivariate}
