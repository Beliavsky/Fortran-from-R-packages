\encoding{latin1}

\name{spatial.sign}
\alias{spatial.sign}

\title{Spatial Signs }
\description{
  Function to obtain the spatial signs of a multivariate dataset. The function can compute the spatial signs also with respect to a given or estimated loacation and scale.
  If both location and scale have to be estimated the \code{\link{HR.Mest}} function is used, if only one has to be estimated the, estimation is done using
  \code{\link{spatial.median}} or \code{\link{tyler.shape}}.
}
\usage{
spatial.sign(X, center = TRUE, shape = TRUE, 
             na.action = na.fail, ...)
}

\arguments{
  \item{X}{ a numeric data frame or matrix.   }
  \item{center}{ either a logical value or a numeric vector of length equal to
          the number of columns of 'X'. See below for more information. }
  \item{shape}{ either a logical value or a square numeric matrix with number of columns equal to
          the number of columns of 'X'. See below for more information.}
  \item{na.action}{ a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.}
  \item{\dots}{ arguments that can be passed on to functions used for the estimation of location
        and shape.}
}
\details{
  The spatial signs U of X with location \eqn{\mu}{mu} and shape V are given by
  \deqn{u_{i}=\frac{(x_{i}-\mu)V^{-\frac{1}{2}}}{\| (x_{i}-\mu)V^{-\frac{1}{2}} \|}.}{u_i = (x_i-mu)V^{-0.5}/||(x_i-mu)V^{-0.5}||.} 
  If a numeric value is given as 'center' and/or 'shape' these are used as \eqn{\mu}{mu} and/or V in the above formula. 
  If 'center' and/or 'shape' are 'TRUE' the values for \eqn{\mu}{mu} and/or V are estimated, if 'FALSE' the origin is used as the 
  value of \eqn{\mu}{mu} and/or the identity matrix as the value of V.

  In the special case of univariate data the univariate signs of the data (centered if requested) are returned and the shape parameter is redundant.
}
\value{
  a matrix with the spatial signs of the data as rows or the univariate signs as a px1 matrix. The
     centering vector and scaling matrix used are returned as attributes
     'center' and 'shape'.

}
\author{Klaus Nordhausen and Seija Sirkia}

\seealso{  \code{\link{HR.Mest}} }
\examples{
set.seed(654321)
cov.matrix <- matrix(c(3,2,1,2,4,-0.5,1,-0.5,2), ncol=3)
X <- rmvnorm(15, c(1,0,-1), cov.matrix)
spatial.sign(X)
spatial.sign(X, center=FALSE, shape=FALSE)
spatial.sign(X, center=colMeans(X), shape=cov(X))
rm(.Random.seed)
}
\keyword{ multivariate }
\keyword{nonparametric}
