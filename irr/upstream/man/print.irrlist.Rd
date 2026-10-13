\name{print.irrlist}
\alias{print.irrlist}
\title{Default printing function for various coefficients of interrater reliability}
\description{
Prints the results of various functions computing coefficients of interrater reliability.
}
\usage{
  \method{print}{irrlist}(x, \dots)
}
\arguments{
  \item{x}{a list with class '"irrlist"' containing the results of the interrater reliability computation.}
  \item{...}{further arguments passed to or from other methods.}
}
\details{
'"print.irrlist"' is only a printing function and is usually not called directly.
}
\author{Matthias Gamer}
\seealso{
  \code{\link[irr]{bhapkar}}, 
  \code{\link[irr]{finn}}, 
  \code{\link[irr]{iota}}, 
  \code{\link[irr]{kappa2}}, 
  \code{\link[irr]{kappam.fleiss}}, 
  \code{\link[irr]{kappam.light}}, 
  \code{\link[irr]{kripp.alpha}}, 
  \code{\link[irr]{kendall}}, 
  \code{\link[irr]{maxwell}}, 
  \code{\link[irr]{meancor}}, 
  \code{\link[irr]{meanrho}}, 
  \code{\link[irr]{rater.bias}}, 
  \code{\link[irr]{robinson}}, 
  \code{\link[irr]{stuart.maxwell}}
}
\examples{
data(anxiety)
# "print.irrlist" is the default printing method of various functions, e.g.
finn(anxiety, 6)
meancor(anxiety)
}
\keyword{print}
