\name{lpochhammer}
\alias{lpochhammer}
\title{ Calculate the logarithm of Pochhammer's symbol }
\description{
  \code{lpochhammer} returns the value of the natural logarithm of Pochhammer's symbol
  calculated as
  
\deqn{\ln \left[ {\left( z \right)_n } \right] = \ln \Gamma \left( {z + n} \right) - \ln \Gamma \left( z \right)}

  where \eqn{ \Gamma \left( z \right) } is the Gamma function
}
\usage{
lpochhammer(z, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{ argument of the symbol }
  \item{n}{ integer number of terms in the symbol }
}
\value{
  The value of the logarithm of the symbol
}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{
  \code{\link{pochhammer}}
}  
\examples{
lpochhammer( pi, 5 )
}
\keyword{ math }
