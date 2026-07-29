\name{ecmDiag}
\alias{ecmDiag}
\title{
Diagnostic Statitics for Symmetric or Asymmetric ECMs
}
\description{
Report a set of diagnostic statistics for symmetric or asymmetric error correction models
}
\usage{
ecmDiag(m, digits = 2)
}
\arguments{
  \item{m}{an object of class \code{ecm} from the function of \code{ecmAsyFit} or \code{ecmSymFit}.}
  \item{digits}{number of digits used in rounding outputs.}
}

\details{
Compute several diagnostic statistics for each ECM equation. This is mainly used to assess the serial correlation in the residuals and model adequacy.
}
\value{
Return a data frame object with the following components by equation: R-squared, Adjusted R-squared, F-statistic, Durbin Watson statistic, p-value for DW statistic, AIC, BIC, and p-value of Ljung_Box Q statistics with 4, 8, 12 autocorrelation coefficients.
}
\references{Enders, W. 2004. Applied Econometric Time Series. John Wiley & Sons, Inc., New York. 480 P.}

\author{Changyou Sun (\email{edwinsun258@gmail.com})}

\seealso{\code{\link{ecmAsyFit}}; \code{\link{ecmSymFit}}; and \code{\link{ecmDiag}}.}

\examples{
# see example at daVich
}
\keyword{regression}