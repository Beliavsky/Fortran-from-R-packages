\name{PcaHubert}
\alias{PcaHubert}
\alias{PcaHubert.formula}
\alias{PcaHubert.default}
\title{ ROBPCA - ROBust method for Principal Components Analysis }
\description{
    The ROBPCA algorithm was proposed by Hubert et al (2005) and stays for
    'ROBust method for Principal Components Analysis'. It is resistant to 
    outliers in the data. The robust loadings are computed using 
    projection-pursuit techniques and the MCD method. Therefore ROBPCA 
    can be applied to both low and high-dimensional data sets. In low 
    dimensions, the MCD method is applied.
}
\usage{
PcaHubert(x, ...)
\method{PcaHubert}{default}(x, k = 0, kmax = 10, alpha = 0.75, mcd = TRUE, skew=FALSE,
maxdir=250, scale = FALSE, signflip = TRUE, crit.pca.distances = 0.975, trace=FALSE, \dots)
\method{PcaHubert}{formula}(formula, data = NULL, subset, na.action, \dots)
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
    of components and then run again. Default is \code{k=0}. }
  \item{kmax}{maximal number of principal components to compute.
    Default is \code{kmax=10}. If \code{k} is provided, \code{kmax} 
    does not need to be specified, unless \code{k} is larger than 10.}
  \item{alpha}{this parameter measures the fraction of outliers the algorithm should
    resist. In MCD alpha controls the size of the subsets over which the 
    determinant is minimized, i.e. alpha*n observations are used for 
    computing the determinant. Allowed values are between 0.5 and 1 
    and the default is 0.75. }
  \item{mcd}{Logical - when the number of variables is sufficiently small,
    the loadings are computed as the eigenvectors of the MCD covariance matrix, 
    hence the function \code{\link{CovMcd}()} is automatically called. The number of 
    principal components is then taken as k = rank(x). Default is \code{mcd=TRUE}.
    If \code{mcd=FALSE}, the ROBPCA algorithm is always applied.}
  \item{skew}{Logical - whether the adjusted outlyingness algorithm for skewed 
    data (Hubert et al., 2009) will be used, default is \code{skew=FALSE}.}
  \item{maxdir}{maximal number of random directions to use for computing the 
    outlyingness (or the adjusted outlyingness when \code{skew=TRUE}) of the data points, 
    see \code{\link[robustbase]{adjOutlyingness}} for more details.. 
    Default is \code{maxdir=250}. If the number \code{n}
    of observations is small all possible \code{n*(n-1)/2} pairs of observations 
    are taken to generate the directions.}
  \item{scale}{a value indicating whether and how the variables should be scaled. 
    If \code{scale=FALSE} (default) or \code{scale=NULL} no scaling is 
    performed (a vector of 1s is returned in the scale slot). If \code{scale=TRUE}
    the data are scaled to have unit variance. Alternatively it can be a function 
    like \code{sd} or \code{mad} or a vector of length equal the number of columns 
    of \code{x}. The value is passed to the underlying function and the result 
    returned is stored in the scale slot. Default is \code{scale=FALSE}.}
  \item{signflip}{a logical value indicating wheather to try to solve the sign indeterminancy of the loadings -   
  ad hoc approach setting the maximum element in a singular vector to be positive. Default is \code{signflip = TRUE}}
  \item{crit.pca.distances}{criterion to use for computing the cutoff 
    values for the orthogonal and score distances. Default is 0.975.}
  \item{trace}{whether to print intermediate results. Default is \code{trace = FALSE}}
}
\details{
\code{PcaHubert}, serving as a constructor for objects of class \code{\link{PcaHubert-class}} 
is a generic function with "formula" and "default" methods. 
The calculation is done using the ROBPCA method of Hubert et al (2005) which can 
be described briefly as follows. For details see the relevant references.

Let \code{n} denote the number of observations, and \code{p} the
number of original variables in the input data matrix \code{X}. The
ROBPCA algorithm finds a robust center \code{M (p x 1)} of the data
and a loading matrix \code{P} which is \code{(p x k)} dimensional.
Its columns are orthogonal and define a new coordinate system. The
scores T, an \code{(n x k)} matrix, are the coordinates of the
centered observations with respect to the loadings:

\deqn{T=(X-M)P} The ROBPCA algorithm also yields a robust covariance
matrix (often singular) which can be computed as

\deqn{S=PLP^t}
where \eqn{L} is the diagonal matrix with the eigenvalues \eqn{l_1, \dots, l_k}.

This is done in the following three main steps:

\bold{Step 1:} The data are preprocessed by reducing their data space to
the subspace spanned by the \code{n} observations. This is done by
singular value decomposition of the input data matrix. As a result
the data are represented using at most \code{n-1=rank(X)} without
loss of information.

\bold{Step 2:} In this step for each data point a measure of outlyingness
is computed. For this purpose the high-dimensional data points are
projected on many univariate directions, each time the univariate
MCD estimator of location and scale is computed  and the
standardized distance to the center is measured. The largest of
these distances (over all considered directions) is the outlyingness
measure of the data point. The \code{h} data points with smallest
outlyingness measure are used to compute the covariance matrix
\eqn{\Sigma_h} and to select the number \code{k} of principal
components to retain. This is done by finding such \code{k} that
\eqn{l_k/l_1 >= 10.E-3} and \eqn{\Sigma_{j=1}^k l_j/\Sigma_{j=1}^r l_j >= 0.8}
Alternatively the number of principal components \code{k} can be
specified by the user after inspecting the scree plot.

\bold{Step 3:} The data points are projected on the k-dimensional subspace
spanned by the \code{k} eigenvectors corresponding to the largest
\code{k} eigenvalues of the matrix \eqn{\Sigma_h}. The location and
scatter of the projected data are computed using the
reweighted MCD estimator and the eigenvectors of this scatter matrix
yield the robust principal components.
}
\value{
  An S4 object of class \code{\link{PcaHubert-class}} which is a subclass of the 
  virtual class \code{\link{PcaRobust-class}}. 
}
\references{
    M. Hubert, P. J. Rousseeuw, K. Vanden Branden (2005), ROBPCA: a new 
    approach to robust principal components analysis, \emph{Technometrics}, \bold{47}, 64--79.   
    
    M. Hubert, P. J. Rousseeuw and T. Verdonck (2009), Robust PCA for skewed data and 
    its outlier map, \emph{Computational Statistics & Data Analysis}, \bold{53}, 2264--2274.

  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} 
}
\note{
    The ROBPCA algorithm is implemented on the bases of the Matlab implementation, 
    available as part of \emph{LIBRA, a Matlab Library for Robust Analysis} to be found at  
    www.wis.kuleuven.ac.be/stat/robust.html
}
\examples{
## PCA of the Hawkins Bradu Kass's Artificial Data
##  using all 4 variables
    data(hbk)
    pca <- PcaHubert(hbk)
    pca

## Compare with the classical PCA
    prcomp(hbk)

## or  
    PcaClassic(hbk)
    
## If you want to print the scores too, use
    print(pca, print.x=TRUE)

## Using the formula interface
    PcaHubert(~., data=hbk)

## To plot the results:

    plot(pca)                    # distance plot
    pca2 <- PcaHubert(hbk, k=2)  
    plot(pca2)                   # PCA diagnostic plot (or outlier map)
    
## Use the standard plots available for prcomp and princomp
    screeplot(pca)    
    biplot(pca)    
    
## Restore the covraiance matrix     
    py <- PcaHubert(hbk)
    cov.1 <- py@loadings \%*\% diag(py@eigenvalues) \%*\% t(py@loadings)
    cov.1    
}
\keyword{robust}
\keyword{multivariate}
