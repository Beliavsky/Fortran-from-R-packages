\name{conc}
\alias{conc}
\alias{Herfindahl}
\alias{Rosenbluth}
\title{Concentration Measures}
\description{
 computes the concentration within a vector according to the
specified concentration measure }
\usage{
conc(x, parameter = NULL, type = c("Herfindahl", "Rosenbluth"), na.rm = TRUE)

Herfindahl(x, parameter = 1, na.rm = TRUE)
Rosenbluth(x, na.rm = TRUE)
}
\arguments{
 \item{x}{a vector containing non-negative elements}
\item{parameter}{parameter of the concentration measure (if set to \code{NULL}
  the default parameter of the respective measure is used)}
\item{type}{character string giving the measure used to compute concentration.
  must be one of the strings in the default argument (the first character
  is sufficient). defaults to "Herfindahl".} 
\item{na.rm}{logical. Should missing values (\code{NA}s) be removed
  prior to computations? If set to \code{FALSE} the computations yield
  \code{NA}.}
}
\details{\code{conc} is just a wrapper for the concentration measures of
         \code{Herfindahl} and \code{Rosenbluth}
         (Hall / Tiedemann / Rosenbluth). If parameter is set to \code{NULL}
         the default from the respective function is used.
}
\value{
  the value of the concentration measure
}
\references{F A Cowell: Measurement of Inequality, 2000, in A B Atkinson / F Bourguignon (Eds): Handbook of Income Distribution, Amsterdam,

F A Cowell: Measuring Inequality, 1995 Prentice Hall/Harvester Wheatshef,

M Hall / N Tidemann: Measures of Concentration, 1967, JASA 62, 162-168.}

\seealso{\code{\link{ineq}}, \code{\link{pov}}}

\examples{
# generate vector (of sales)
x <- c(541, 1463, 2445, 3438, 4437, 5401, 6392, 8304, 11904, 22261)
# compute Herfindahl coefficient with parameter 1
conc(x)
# compute coefficient of Hall/Tiedemann/Rosenbluth
conc(x, type="Rosenbluth")
}
\keyword{univar}

