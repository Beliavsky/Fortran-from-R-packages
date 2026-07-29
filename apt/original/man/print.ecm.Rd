\name{print.ecm}
\alias{print.ecm}
\title{Printing Outputs from Error Correction Models}
\description{
Show main outputs from symmetric \code{ecmSymFit} or asymmetric \code{ecmAsyFit} error correction models.
}
\usage{\method{print}{ecm}(x, \dots)}
\arguments{
  \item{x}{an object of class \code{ecm} from the function of \code{ecmAsyFit} or \code{ecmSymFit}.}
  \item{\dots}{additional arguments to be passed.}
}
\details{This is the print method for \code{ecmAsyFit} or \code{ecmSymFit} to show main model outputs.}
\value{
Summary results of the two ECM equations are shown for the two focal price series.
}
\author{Changyou Sun (\email{edwinsun258@gmail.com})}

\seealso{\code{\link{ecmSymFit}} and \code{\link{ecmAsyFit}}.}

\examples{
# see example at daVich
}
\keyword{methods}
