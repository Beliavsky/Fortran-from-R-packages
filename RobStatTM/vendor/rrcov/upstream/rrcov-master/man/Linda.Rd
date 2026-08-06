\name{Linda}
\alias{Linda}
\alias{Linda.formula}
\alias{Linda.default}
\title{ Robust Linear Discriminant Analysis }
\description{
    Robust linear discriminant analysis based on MCD and returns
    the results as an object of class \code{Linda} (aka constructor).
}
\usage{
Linda(x, ...)

\method{Linda}{default}(x, grouping, prior = proportions, tol = 1.0e-4,
                 method = c("mcd", "mcdA", "mcdB", "mcdC", "fsa", "mrcd", "ogk"),
                 alpha=0.5, l1med=FALSE, cov.control, trace=FALSE, ...)
}
\arguments{
  \item{x}{a matrix or data frame containing the explanatory variables (training set). }
  \item{grouping}{grouping variable:  a factor specifying the class for each observation.}
  \item{prior}{prior probabilities, default to the class proportions for the training set.}
  \item{tol}{tolerance}
  \item{method}{method}
  \item{alpha}{this parameter measures the fraction of outliers the algorithm should
    resist. In MCD alpha controls the size of the subsets over which the
    determinant is minimized, i.e. alpha*n observations are used for
    computing the determinant. Allowed values are between 0.5 and 1
    and the default is 0.5.}
  \item{l1med}{whether to use L1 median (space median) instead of MCD to compute
    the group means locations in order to center the data in methods \code{mcdB} and \code{mcdC}. 
    useful in case of groups with small size. Default is \code{l1med = FALSE}.}
  \item{cov.control}{specifies which covariance estimator to use by providing 
    a \code{\link{CovControl-class}} object. 
    The default is \code{\link{CovControlMcd-class}} which will 
    indirectly call \code{\link{CovMcd}}. 
    If \code{cov.control=NULL} is specified, the classical estimates 
    will be used by calling \code{\link{CovClassic}}}.
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}.}
  \item{\dots}{arguments passed to or from other methods}
}
\details{
  details
}
\value{
  Returns an S4 object of class \code{Linda}
}
\references{
    Hawkins, D.M. and McLachlan, G.J. (1997) High-Breakdown Linear Discriminant Analysis,
    \emph{Journal of the American Statistical Association}, \bold{92}, 136--143.

    Todorov V. (2007) Robust selection of variables in linear discriminant analysis,
        \emph{Statistical Methods and Applications}, \bold{15}, 395--407, doi:10.1007/s10260-006-0032-6.

    Todorov, V. and Pires, A.M. (2007) Comparative Performance of Several Robust
        Linear Discriminant Analysis Methods.
        \emph{REVSTAT Statistical Journal}, \bold{5}, p 63--83.

  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at}
}
\seealso{  \code{\link[rrcov]{CovMcd}}, \code{\link[rrcov]{CovMrcd}}  }
\examples{
## Example anorexia
library(MASS)
data(anorexia)

## start with the classical estimates
lda <- LdaClassic(Treat~., data=anorexia)
predict(lda)@classification

## try now the robust LDA with the default method (MCD with pooled whitin cov matrix)
rlda <- Linda(Treat~., data= anorexia)
predict(rlda)@classification

## try the other methods
Linda(Treat~., data= anorexia, method="mcdA")
Linda(Treat~., data= anorexia, method="mcdB")
Linda(Treat~., data= anorexia, method="mcdC")

## try the Hawkins&McLachlan method
## use the default method
grp <- anorexia[,1]
grp <- as.factor(grp)
x <- anorexia[,2:3]
Linda(x, grp, method="fsa")

## Do DA with Linda and method mcdB or mcdC, when some classes
## have very few observations. Use L1 median instead of MCD
##  to compute the group means (l1med=TRUE).

data(fish)

# remove observation #14 containing missing value
fish <- fish[-14,]

# The height and width are calculated as percentages 
#   of the third length variable
fish[,5] <- fish[,5]*fish[,4]/100
fish[,6] <- fish[,6]*fish[,4]/100

table(fish$Species) 
Linda(Species~., data=fish, l1med=TRUE)
Linda(Species~., data=fish, method="mcdC", l1med=TRUE)

}
\keyword{robust}
\keyword{multivariate}
