\name{summary.armasubsets}
\alias{summary.armasubsets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary of output from the armasubsets function}
\description{
Add the calculation of AIC and AICc. See the help manual of
regsubsets function of the leaps package}
\usage{
\method{summary}{armasubsets}(object, all.best = TRUE, matrix = TRUE, matrix.logical = FALSE,
 df = NULL, ...)
}
\arguments{
\item{object}{armasubsets object}
\item{all.best}{Show all the best subsets or just one of each size}
\item{matrix}{Show a matrix of the variables in each model or just summary statistics}
\item{matrix.logical}{With matrix=TRUE, the matrix is logical TRUE/FALSE or string "*"/code{" "}}
\item{df}{Specify a number of degrees of freedom for the summary statistics. The default is n-1}
\item{...}{Other arguments for future methods}
}
\author{ Kung-Sik Chan, based on previous work of Thomas Lumley}
\keyword{methods}
