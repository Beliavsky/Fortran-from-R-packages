\name{fmin.gridsearch}
\alias{fmin.gridsearch}
\title{
  Grid evaluation of an unconstrained cost function
}
\description{
  Evaluate an unconstrained cost function on a grid of points around a given 
  initial point estimate.
}
\usage{
  fmin.gridsearch(fun = NULL, x0 = NULL, xmin = NULL, 
                  xmax = NULL, npts = 3, alpha = 10)
}
\arguments{
  \item{fun}{An unconstrained cost function returning a numeric scalar, similar 
    to those used in the \code{fminsearch} function.}
  \item{x0}{The initial point estimate, provided as a numeric vector.}
  \item{xmin}{Optional: a vector of lower bounds.}
  \item{xmax}{Optional: a vector of upper bounds.}
  \item{npts}{An integer scalar greater than 2, indicating the number of 
    evaluation points will be used on each dimension to build the search grid.}
  \item{alpha}{A vector of numbers greater than 1, which give the factor(s) used
    to calculate the evaluation range of each dimension of the search grid (see 
    Details). If \code{alpha} length is lower than that of \code{x0}, elements 
    of \code{alpha} are recycled. If its length is higher than that of 
    \code{x0}, \code{alpha} is truncated.}
}
\details{
  \code{fmin.gridsearch} evaluates the cost function at each point 
  of a grid of \code{npts^length(x0)} points. If lower (\code{xmin}) and upper 
  (\code{xmax}) bounds are provided, the range of evaluation points is limited 
  by those bounds and \code{alpha} is not used. Otherwise, the range of 
  evaluation points is defined as \code{[x0/alpha,x0*alpha]}.
  
  The actual evaluation of the cost function is delegated to 
  \code{optimbase.gridsearch}.
}
\value{
  Return a data.frame with the coordinates of the evaluation point, the value of
  the cost function and its feasibility. Because the cost function is 
  unconstrained, it is always feasible. The data.frame is ordered by feasibility
  and increasing value of the cost function.
}
\author{
  Sebastien Bihorel (\email{sb.pmlab@gmail.com})
}
\seealso{
  \code{\link[neldermead]{fminsearch}}, 
  \code{\link[optimbase]{optimbase.gridsearch}}
}
\keyword{ method }