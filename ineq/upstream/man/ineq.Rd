\name{ineq}
\alias{ineq}
\alias{Gini}
\alias{RS}
\alias{Atkinson}
\alias{Theil}
\alias{Kolm}
\alias{var.coeff}
\alias{entropy}
\title{Inequality Measures}
\description{
 computes the inequality within a vector according to the
specified inequality measure }
\usage{
ineq(x, parameter = NULL, type = c("Gini", "RS", "Atkinson", "Theil", "Kolm", "var",
     "square.var", "entropy"), na.rm = TRUE)

Gini(x, corr = FALSE, na.rm = TRUE)
RS(x, na.rm = TRUE)
Atkinson(x, parameter = 0.5, na.rm = TRUE)
Theil(x, parameter = 0, na.rm = TRUE)
Kolm(x, parameter = 1, na.rm = TRUE)
var.coeff(x, square = FALSE, na.rm = TRUE)
entropy(x, parameter = 0.5, na.rm = TRUE)
}
\arguments{
 \item{x}{a vector containing at least non-negative elements}
\item{parameter}{parameter of the inequality measure (if set to \code{NULL}
  the default parameter of the respective measure is used)}
\item{type}{character string giving the measure used to compute inequality.
  must be one of the strings in the default argument (the first character
  is sufficient). defaults to "Gini".}
\item{corr}{logical. Argument of the function \code{Gini} specifying whether
  or not a finite sample correction should be applied.}
\item{square}{logical. Argument of the function \code{var.coeff}, for details
  see below.}
\item{na.rm}{logical. Should missing values (\code{NA}s) be removed
  prior to computations? If set to \code{FALSE} the computations yield
  \code{NA}.}
}
\details{\code{ineq} is just a wrapper for the inequality measures \code{Gini},
         \code{RS}, \code{Atkinson}, \code{Theil}, \code{Kolm},\code{var.coeff},
         \code{entropy}. If parameter is set to \code{NULL} the default from
         the respective function is used.

         \code{Gini} is the Gini coefficient, \code{RS} is the the Ricci-Schutz
         coefficient (also called Pietra's measure), \code{Atkinson} gives
         Atkinson's measure and \code{Kolm} computes Kolm's measure.

         If the parameter in \code{Theil} is 0 Theil's entropy measure is
         computed, for every other value Theil's second measure is
         computed.

         \code{ineq(x, type="var")} and \code{var.coeff(x)} respectively
         compute the coefficient of variation, while
         \code{ineq(x,type="square.var")} and \code{var.coeff(x, square=TRUE)}
         compute the squared coefficient of variation.

         \code{entropy} computes the generalized entropy, which is for
         parameter 1 equal to Theil's entropy coefficient and for parameter
         0 equal to the second measure of Theil.
}
\value{
  the value of the inequality measure
}
\references{F A Cowell: Measurement of Inequality, 2000, in A B Atkinson / F Bourguignon (Eds): Handbook of Income Distribution, Amsterdam,

F A Cowell: Measuring Inequality, 1995 Prentice Hall/Harvester Wheatshef,

Marshall / Olkin: Inequalities: Theory of Majorization and Its
Applications, New York 1979 (Academic Press).}

\seealso{\code{\link{conc}}, \code{\link{pov}}}

\examples{
# generate vector (of incomes)
x <- c(541, 1463, 2445, 3438, 4437, 5401, 6392, 8304, 11904, 22261)
# compute Gini coefficient
ineq(x)
# compute Atkinson coefficient with parameter=0.5
ineq(x, parameter=0.5, type="Atkinson")
}
\keyword{univar}

