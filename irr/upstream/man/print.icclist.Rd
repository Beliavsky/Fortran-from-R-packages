\name{print.icclist}
\alias{print.icclist}
\title{Default printing function for ICC results}
\description{
Prints the results of the ICC computation.
}
\usage{
  \method{print}{icclist}(x, \dots)
}
\arguments{
  \item{x}{a list with class '"icclist"' containing the results of the ICC computation.}
  \item{...}{further arguments passed to or from other methods.}
}
\details{
'"print.icclist"' is only a printing function and is usually not called directly.
}
\author{Matthias Gamer}
\seealso{
  \code{\link[irr]{icc}}
}
\examples{
data(anxiety)
# "print.icclist" is the default printing function of "icc"
icc(anxiety, model="twoway", type="agreement")
}
\keyword{print}
