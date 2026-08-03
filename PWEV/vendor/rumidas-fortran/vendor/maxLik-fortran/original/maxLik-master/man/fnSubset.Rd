\name{fnSubset}
\alias{fnSubset}
\title{
  Call fnFull with variable and fixed parameters
}
\description{
  Combine variable parameters with with fixed parameters and pass to
  \code{fnFull}.  Useful for optimizing over a subset of parameters
  without writing a separate function.  Values are combined by name if
  available.  Otherwise, \code{xFull} is constructed by position (the
  default).  
}
\usage{
fnSubset(x, fnFull, xFixed, xFull=c(x, xFixed), ...) 
}
\arguments{
  \item{x}{
    Variable parameters to be passed to \code{fnFull}.  
  }
  \item{fnFull}{
    Function whose first argument has length = length(xFull).  
  }
  \item{xFixed}{
    Parameter values to be combined with \code{x} to construct the first
    argument for a call to \code{fnFull}.  
  }
  \item{xFull}{
    Prototype initial argument for \code{fnFull}.  
  }
  \item{\dots}{
    Optional arguments passed to \code{fnFull}.  
  }
}
\details{
  This function first confirms that
  \code{length(x) + length(xFixed) == length(xFull)}.
  Next,
  \itemize{
    \item If \code{xFull} has names, match at least \code{xFixed} by
    name.
    \item Else \code{xFull = c(x, xFixes)}, the default.  
  }
  Finally, call \code{fnFull(xFull, ...)}.
}  
\value{
  value returned by \code{fnFull} 
}
%\references{ }
\author{
  Spencer Graves 
}
\seealso{
  \code{\link{optim}}
  \code{\link[dlm]{dlmMLE}}
  \code{\link{maxLik}}
  \code{\link{maxNR}} 
}
\examples{
##
## Example with 'optim'
##
fn <- function(x) (x[2]-2*x[1])^2
# note: true minimum is 0 on line 2*x[1] == x[2]
fullEst <- optim(par=c(1,1), method="BFGS", fn=fn)
fullEst$par
# par = c(0.6, 1.2) at minimum (not convex)

# Fix the last component to 4 
est4 <- optim(par=1, fn=fnSubset, method="BFGS", fnFull=fn, xFixed=4)
est4$par
# now there is a unique minimun x[1] = 2

# Fix the first component
fnSubset(x=1, fnFull=fn, xFixed=c(a=4), xFull=c(a=1, b=2))
# After substitution:  xFull = c(a=4, b=1),
# so fn = (1 - 2*4)^2 = (-7)^2 = 49

est4. <- optim(par=1, fn=fnSubset, method="BFGS",
               fnFull=fn, xFixed=c(a=4), 
               xFull=c(a=1, b=2))
est4.$par
# At optimum: xFull=c(a=4, b=8),
# so fn = (8 - 2*4)^2 = 0

##
## Example with 'maxLik'
##
fn2max <- function(x) -(x[2]-2*x[1])^2
# -> need to have a maximum
max4 <- maxLik(fnSubset, start=1, fnFull=fn2max, xFixed=4)
summary(max4)
# Similar result using fixed parameters in maxNR, called by maxLik 
max4. <- maxLik(fn2max, start=c(1, 4), fixed=2)
summary(max4.)
}
\keyword{optimize}
\keyword{utilities}
