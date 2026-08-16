\name{pochhammer}
\alias{pochhammer}
\title{ Calculate the value of Pochhammer's symbol }
\description{
  \code{pochhammer} returns the value of Pochhammer's symbol calculated as
  
  \deqn{\left( z \right)_n  = z\;\left( {z + 1} \right)\; \ldots \;\left( {z + n - 1} \right) = \frac{{\Gamma \left( {z + n} \right)}}{{\Gamma \left( z \right)}}}
  
  where \eqn{\Gamma \left( z \right)} is the Gamma function

}
\usage{
pochhammer(z, n)
}
\arguments{
  \item{z}{ numeric value for the argument of the symbol }
  \item{n}{ integer value for the number of terms in the symbol }
}
\value{
  The value of Pochhammer's symbol
}
\author{ Frederick Novomestky \email{fnovomes@poly.edu}
}
\examples{
###
### compute the Pochhamer's symbol fo z equal to 1 and
### n equal to 5
###
pochhammer( 1, 5 )
}
\keyword{ math }
