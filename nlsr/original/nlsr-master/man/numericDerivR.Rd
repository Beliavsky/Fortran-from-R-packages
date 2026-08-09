\name{numericDerivR}
\alias{numericDerivR}
\title{numericDerivR: numerically evaluates the gradient of an expression. All in R}
\usage{
numericDerivR(
  expr,
  theta,
  rho = parent.frame(),
  dir = 1,
  eps = .Machine$double.eps^(1/if (central) 3 else 2),
  central = FALSE
)
}
\arguments{
\item{expr}{expression or call to be differentiated. Should evaluate to a numeric vector.}

\item{theta}{character vector of names of numeric variables used in expr.}

\item{rho}{environment containing all the variables needed to evaluate expr.}

\item{dir}{numeric vector of directions, typically with values in -1, 1 to use for the finite 
differences; will be recycled to the length of theta.}

\item{eps}{a positive number, to be used as unit step size hh for the approximate 
numerical derivative (f(x+h)-f(x))/h (f(x+h)-f(x))/h or the central version, see central.}

\item{central}{logical indicating if central divided differences should be computed, 
i.e., (f(x+h) - f(x-h)) / 2h (f(x+h)-f(x-h))/2h. These are typically more accurate but 
need more evaluations of f()f().}
}
\value{
The value of eval(expr, envir = rho) plus a matrix attribute "gradient". The columns of 
this matrix are the derivatives of the value with respect to the variables listed in theta.
}
\description{
This version is all in R to replace the C version in package stats
}
\examples{
ex <- expression(a/(1+b*exp(-tt*c)) - weed)
weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
              38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12
a <- 200; b <- 50; c <- 0.3
dhobb <- numericDerivR(ex, theta=c("a", "b", "c"))
print(dhobb)
# exf <- ~ a/(1+b*exp(-tt*c)) - weed
# Note that a formula doesn't work
# dh1 <-  try(numericDerivR(exf, theta=c("a", "b", "c")))
}
