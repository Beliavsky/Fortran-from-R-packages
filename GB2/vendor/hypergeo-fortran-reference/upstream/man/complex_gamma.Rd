\name{complex_gamma}
\alias{complex_gamma}
\alias{complex_factorial}
\alias{lanczos}
\title{
Gamma function for complex arguments
}
\description{
Gamma and factorial functions for complex arguments
}
\usage{
complex_gamma(z, log = FALSE)
complex_factorial(z, log = FALSE)
lanczos(z,log=FALSE)
}
\arguments{
  \item{z}{Primary argument, a complex vector}
  \item{log}{Boolean, with default \code{FALSE} meaning to return the
    function value and \code{TRUE} meaning to return its logarithm}
}
\details{
 Method follows that of Lanczos, coefficients identical to those of the GSL
}
\references{
Lanczos, C. 1964. \dQuote{A precision approximation of the gamma
  function}. \emph{Journal of the society for industrial and applied
  mathematics series B}, Volume 1, pp86-96

 M. Galassi et al, GNU Scientific Library Reference Manual (3rd Ed.), ISBN 0954612078. 

}
\author{
Robin K. S. Hankin
}

\examples{

complex_gamma(5)  # should be 4!=24

complex_gamma(1+1i)   # takes complex arguments

complex_gamma(-5/2) + sqrt(pi)*8/15  # should be small

z <- pi + 1i*sqrt(2)
complex_gamma(z+1)-z*complex_gamma(z)   # should be small


complex_gamma(z)*complex_gamma(1-z) - pi/sin(pi*z)  # small


}