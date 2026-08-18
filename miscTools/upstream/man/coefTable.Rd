\name{coefTable}
\alias{coefTable}
\title{Coefficient Table}

\description{
   Generate Table for Coefficients, Std. Errors, t-values and P-values.
}

\usage{
   coefTable( coef, stdErr, df = NULL )
}

\arguments{
   \item{coef}{vector that contains the coefficients.}
   \item{stdErr}{vector that contains the standard errors of the coefficients.}
   \item{df}{degrees of freedom of the t-test used to calculate P-values.}
}

\value{
   a matrix with 4 columns: coefficients, standard errors, t-values
      and P-values.
   If argument \code{df} is not provided, the last column (P-values)
   is filled with \code{NA}s.
}

\author{Arne Henningsen}

\examples{
   coefTable( rnorm( 10 ), 0.5 * abs( rnorm( 10 ) ), 20 )
}

\keyword{models}
