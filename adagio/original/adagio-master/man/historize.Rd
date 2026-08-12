\name{Historize}
\alias{Historize}
\title{
  Historize function
}
\description{
  Provides storage for function calls.
}
\usage{
Historize(fun, len = 0, ...)
}
\arguments{
  \item{fun}{Function of one or several variables.}
  \item{len}{If > 0, input values will be stored, too.}
  \item{...}{Additional parameters to be handed to the function.}
}
\details{
  Historize() adds storage to the function. If \code{len=0} then only
  function values will be stored in a vector, if \code{len>0} then
  variables will be stored in a vector, and function values will be stored
  in a matrix, one line for every function call.

  If \code{Fn = Historize(...)} is the 'historized' function, then the
  storage can be retrieved with \code{Fn()}, and can be cleared with
  \code{Fn(NULL)}.

  Filling the storage will take extra time and can slow down the
  function calls. Especially also storing the variables used in the
  call (with \code{len>0}) will make it considerably slower.

  Functions can have multivariate output; the user is asked to take
  care of handling the output vector or matrix correctly. The function
  may even require additional parameters.
}
\value{
  Returns a list, $input the input variables as matrix, $H2 the function
  values as vector, $nvars the number of input variables of the function,
  and $ncalls the number or recorded function calls.
}
\author{
  Hans W. Borchers
}
\note{
  Can also be applied to functions that output a vector (same length
  for every call).
}
\seealso{
  \code{trace}
}
\examples{
f <- function(x) sum(x^2)
F <- Historize(f, len = 1)
c( F(c(1,1)), F(c(1,2)), F(c(2,1)), F(c(2,2)) )
#> [1] 2 5 5 8
F()
#> $input
#>      [,1] [,2]
#> [1,]    1    1
#> [2,]    1    2
#> [3,]    2    1
#> [4,]    2    2
#> 
#> $values
#> [1] 2 5 5 8
#> 
#> $nvars
#> [1] 2
#> 
#> $ncalls
#> [1] 4
F(NULL)     # reset memory


## Rastrigin under Differential Evolution

Fn <- Historize(fnRastrigin)
dm <- 10
lb <- rep(-5.2, dm); ub <- -lb

sol <- simpleDE(Fn, lower = lb, upper = ub)
fvalues <- Fn()$values
fvals <- cummin(fvalues)

\dontrun{
plot(fvalues, type = 'p', col = 7, pch = '.', cex = 2,
     main = "Simple DE optimization", xlab = '', ylab = '')
lines(fvals, col = 2, lwd = 2)
legend(length(fvalues), max(fvalues),
       c("Intermediate values", "cummulated min"),
       xjust = 1, col = c(7, 2), lwd = 2)
grid()}
}
\keyword{ manip }
