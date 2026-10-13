\name{pov}
\alias{pov}
\alias{Watts}
\alias{Sen}
\alias{SST}
\alias{Foster}
\title{Poverty Measures}
\description{
 computes the poverty of an (income) vector according to the
specified poverty measure }
\usage{
pov(x, k, parameter = NULL, type = c("Watts", "Sen", "SST", "Foster"), na.rm = TRUE)

Watts(x, k, na.rm = TRUE)
Sen(x, k, na.rm = TRUE)
SST(x, k, na.rm = TRUE)
Foster(x, k, parameter = 1, na.rm = TRUE)
}
\arguments{
 \item{x}{a vector containing at least non-negative elements}
 \item{k}{a constant giving the absolute poverty line}
 \item{parameter}{parameter of the poverty measure (if set to \code{NULL}
  the default parameter of the respective measure is used)}
 \item{type}{character string giving the measure used to compute poverty
    coefficient must be one of the strings in the default argument. Defaults to "Watts".}
 \item{na.rm}{logical. Should missing values (\code{NA}s) be removed
   prior to computations? If set to \code{FALSE} the computations yield
   \code{NA}.}
}
\details{\code{pov} is just a wrapper for the poverty measures of
         \code{Watts}, \code{Sen}, \code{SST}, and \code{Foster}
         (Foster / Greer / Thorbecke). If parameter is set to \code{NULL} the
         default from the respective function is used.

         \code{Foster} gives for parameter 1 the headcount ratio and for
         parameter 2 the poverty gap ratio.
}
\value{
  the value of the poverty measure
}
\references{
  Foster, J. E. (1984). On Economic Poverty: A Survey of Aggregate Measures.
    \emph{Advances in Econometrics}, \bold{3}, 215--251. 

  Shorroks, A. F. (1995). Revisiting the Sen Poverty Index.
    \emph{Econometrica}, \bold{63}(5), 1225--1230.

  Zheng, B. (1997). Aggregate Poverty Measures.
    \emph{Journal of Economic Surveys}, \bold{11}, 123--162.  
}

\seealso{\code{\link{ineq}}, \code{\link{conc}}}

\examples{
# generate vectors (of incomes)
x <- c(541, 1463, 2445, 3438, 4437, 5401, 6392, 8304, 11904, 22261)
y <- c(841, 2063, 2445, 3438, 4437, 5401, 6392, 8304, 11304, 21961)
# compute Watts index with poverty line 2000
pov(x, 2000)
pov(y, 2000)
# compute headcount ratio with poverty line 2000
pov(x, 2000, parameter=1, type="Foster")
pov(y, 2000, parameter=1, type="Foster")}
\keyword{univar}
