\name{summary.maxim}
\alias{summary.maxim}
\alias{print.summary.maxim}
\title{Summary method for maximization}
\description{
  Summarizes the general maximization results in a way that does not
  assume the function is log-likelihood.
}
\usage{
\method{summary}{maxim}( object, hessian=FALSE, unsucc.step=FALSE, ... )
\method{print}{summary.maxim}(x,
                              max.rows=getOption("max.rows", 20),
                              max.cols=getOption("max.cols", 7),
                              ... )
}
\arguments{
  \item{object}{optimization result, object of class
    \code{maxim}.  See \code{\link{maxNR}}.}
  \item{hessian}{logical, whether to display Hessian matrix.}
  \item{unsucc.step}{logical, whether to describe last unsuccesful step
    if \code{code} == 3}
  \item{x}{object of class \code{summary.maxim}, summary of maximization
    result.
  }
  \item{max.rows}{maximum number of rows to be printed.  This applies to
    the resulting coefficients (as those are printed as a matrix where
    the other column is the gradient), and to the Hessian if requested.
  }
  \item{max.cols}{maximum number of columns to be printed.  Only Hessian
    output, if requested, uses this argument.
  }
  \item{\ldots}{currently not used.}
}
\value{
  Object of class \code{summary.maxim}, intended to be printed with
  corresponding print method.
}
\author{Ott Toomet}
\seealso{\code{\link{maxNR}}, \code{\link{returnCode}},
  \code{\link{returnMessage}}} 
\examples{
## minimize a 2D quadratic function:
f <- function(b) {
  x <- b[1]; y <- b[2];
  val <- -(x - 2)^2 - (y - 3)^2  # concave parabola
  attr(val, "gradient") <- c(-2*x + 4, -2*y + 6)
  attr(val, "hessian") <- matrix(c(-2, 0, 0, -2), 2, 2)
  val
}
## Note that NR finds the minimum of a quadratic function with a single
## iteration.  Use c(0,0) as initial value.  
res <- maxNR( f, start = c(0,0) ) 
summary(res)
summary(res, hessian=TRUE)
}
\keyword{methods}
\keyword{print}
