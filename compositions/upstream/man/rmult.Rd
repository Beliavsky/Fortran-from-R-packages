\name{rmult}
\alias{rmult}
\alias{print.rmult}
\title{Simple treatment of real vectors }
\description{
  A class to collect real multivariate vectors.
}
\usage{
 rmult(X,parts=1:NCOL(oneOrDataset(X)),orig=gsi.orig(X),
                missingProjector=attr(X,"missingProjector"),
                V = gsi.getV(X))
 \method{print}{rmult}(x,..., verbose=FALSE)
          }
\arguments{
  \item{X}{vector or dataset of numbers considered as elements of a R-vector}
  \item{parts}{vector containing the indices  xor names of the columns
    to be used}
  \item{x}{an rmult object}
  \item{orig}{the original untransformed dataset}
  \item{missingProjector}{the Projector on the observed subspace}
  \item{V}{the \emph{inverse} of the transformation matrix}
  \item{\dots}{further generic arguments passed to \code{print.default}}
  \item{verbose}{logical, do you want to get all information about original data and transformation function (if any) with a \code{print} call? defaults to FALSE (to print strict content only)}
}
\value{
  a vector of class \code{"rmult"} representing one vector
  or a matrix of class \code{"rmult"}, representing
  multiple vectors by rows.  
}
\details{
The \code{rmult} class is a simple convenience class to treat
data in the scale of real vectors just like data in the scale of real
numbers. A major aspect to take into account is that the internal arithmetic of R is
different for these vectors, e.g. \code{mean} works as \code{colMeans} in a data frame,
or matrix-vector operations are done row-wise.    
}

\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{+.rmult}}, \code{\link{scalar}}, \code{\link{norm.rmult}},
  \code{\link{\%*\%.rmult}},
  \code{\link{rplus}}, \code{\link{acomp}},
}

\examples{
plot(rnorm.rmult(30,mean=0:4,var=diag(1:5)+10))

}
\keyword{multivariate}
