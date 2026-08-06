\name{Pca-class}
\docType{class}
\alias{Pca-class}
\alias{getCenter,Pca-method}
\alias{getScale,Pca-method}
\alias{getEigenvalues,Pca-method}
\alias{getLoadings,Pca-method}
\alias{getPrcomp,Pca-method}
\alias{getScores,Pca-method}
\alias{getSdev,Pca-method}
\alias{plot,Pca,missing-method}
\alias{show,Pca-method}
\alias{predict,Pca-method}
\alias{predict}
\alias{screeplot,Pca-method}
\alias{screeplot}
\alias{summary,Pca-method}

\title{Class "Pca" - virtual base class for all classic and robust PCA classes 
}
\description{  The class \code{Pca} searves as a base class for deriving all other 
classes representing the results of the classical and robust Principal 
Component Analisys methods
}
\section{Objects from the Class}{A virtual Class: No objects may be created from it.}
\section{Slots}{
     \describe{
    \item{\code{call}:}{Object of class \code{"language"} }
    \item{\code{center}:}{Object of class \code{"vector"} the center of the data }
    \item{\code{scale}:}{Object of class \code{"vector"} the scaling applied to each variable of the data }
    \item{\code{rank}:}{Object of class \code{"numeric"} the rank of the data matrix }
    \item{\code{loadings}:}{Object of class \code{"matrix"} the matrix 
    of variable loadings (i.e., a matrix whose columns contain the eigenvectors) }
    \item{\code{eigenvalues}:}{Object of class \code{"vector"} the eigenvalues }
    \item{\code{scores}:}{Object of class \code{"matrix"} the scores - the value 
    of the projected on the space of the principal components data (the centred 
    (and scaled if requested) data multiplied 
    by the \code{loadings} matrix) is returned.  Hence, \code{cov(scores)} 
    is the diagonal matrix \code{diag(eigenvalues)} }
    \item{\code{k}:}{Object of class \code{"numeric"} number of (choosen) principal components }
    \item{\code{sd}:}{Object of class \code{"Uvector"} Score distances within the robust PCA subspace }
    \item{\code{od}:}{Object of class \code{"Uvector"} Orthogonal distances to the robust PCA subspace }
    \item{\code{cutoff.sd}:}{Object of class \code{"numeric"} Cutoff value for the score distances}
    \item{\code{cutoff.od}:}{Object of class \code{"numeric"} Cutoff values for the orthogonal distances }
    \item{\code{flag}:}{Object of class \code{"Uvector"} The observations whose score distance is larger 
    than cutoff.sd or whose orthogonal distance is larger than cutoff.od can be considered 
    as outliers and receive a flag equal to zero.
    The regular observations receive a flag 1 }
    \item{crit.pca.distances}{criterion to use for computing the cutoff values for the orthogonal and score distances. Default is 0.975.}
    \item{\code{n.obs}:}{Object of class \code{"numeric"} the number of observations }
    \item{\code{eig0}:}{Object of class \code{"vector"} all eigenvalues }
    \item{\code{totvar0}:}{Object of class \code{"numeric"} the total variance explained (\code{=sum(eig0)}) }
  }
}
\section{Methods}{
  \describe{
    \item{getCenter}{\code{signature(obj = "Pca")}: center of the data }
    \item{getScale}{\code{signature(obj = "Pca")}: return the scaling applied to each variable }
    \item{getEigenvalues}{\code{signature(obj = "Pca")}: the eigenvalues of the 
        covariance/correlation matrix, though the calculation is actually done 
        with the singular values of the data matrix) }
    \item{getLoadings}{\code{signature(obj = "Pca")}: returns the matrix 
        \code{loadings} (i.e., a matrix whose columns contain the eigenvectors). 
        The function prcomp returns this matrix in the element rotation. }
    \item{getPrcomp}{\code{signature(obj = "Pca")}: returns an S3 object \code{prcomp} 
        for compatibility with the functions prcomp() and princomp(). Thus the 
        standard plots screeplot() and biplot() can be used}
    \item{getScores}{\code{signature(obj = "Pca")}:  returns the rotated data (the centred 
        (and scaled if requested) data multiplied by the loadings matrix). }
    \item{getSdev}{\code{signature(obj = "Pca")}: returns the standard deviations of the 
        principal components (i.e., the square roots of the eigenvalues of the 
        covariance/correlation matrix, though the calculation is actually done 
        with the singular values of the data matrix) }
    \item{plot}{\code{signature(x = "Pca")}: produces a distance plot (if \code{k=rank}) or 
        distance-distance plot (if\code{k<rank}) }
    \item{print}{\code{signature(x = "Pca")}: prints the results. The difference to the show()
        method is that additional parametesr are possible.}
    \item{show}{\code{signature(object = "Pca")}: prints the results }
    \item{predict}{\code{signature(object = "Pca")}: calculates prediction using the results in 
        \code{object}. An optional data frame or matrix in which to look for variables with which 
        to predict. If omitted, the scores are used. If the original fit used a formula or 
        a data frame or a matrix with column names, newdata must contain columns with the 
        same names. Otherwise it must contain the same number of columns, 
        to be used in the same order.  See also \code{\link{predict.prcomp}} and 
        \code{\link{predict.princomp}}}
    \item{screeplot}{\code{signature(x = "Pca")}: plots the variances against the 
        number of the principal component. See also \code{\link{plot.prcomp}} and 
        \code{\link{plot.princomp}}}
  }
}
\references{
  Todorov V & Filzmoser P (2009), An Object Oriented Framework for Robust Multivariate Analysis. 
  \emph{Journal of Statistical Software}, \bold{32}(3), 1--47.
  \doi{10.18637/jss.v032.i03}.
}
\author{ Valentin Todorov \email{valentin.todorov@chello.at} }
\seealso{
     \code{\link{PcaClassic}}, \code{\link{PcaClassic-class}}, \code{\link{PcaRobust-class}}
}
\examples{
showClass("Pca")
}
\keyword{robust}
\keyword{multivariate}
