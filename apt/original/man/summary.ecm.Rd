\name{summary.ecm}
\alias{summary.ecm}
\title{Summary of Results from Error Correction Model}
\description{This summarizes the main results from error correction models.}
\usage{\method{summary}{ecm}(object, digits=3, \dots)}
\arguments{
  \item{object}{an object of class \code{ecm} from the function of \code{ecmAsyFit} or \code{ecmSymFit}.}
  \item{digits}{number of digits for rounding outputs}
  \item{\dots}{additional arguments to be passed.}
}
\details{This wraps up the coefficents and statistics from ECM by equation.}
\value{A data frame object with coefficients and related statistics by equation.}
\author{Changyou Sun (\email{edwinsun258@gmail.com})}

\seealso{\code{\link{ecmSymFit}} and \code{\link{ecmAsyFit}}.}
\examples{
# see example at daVich
}
\keyword{methods}
