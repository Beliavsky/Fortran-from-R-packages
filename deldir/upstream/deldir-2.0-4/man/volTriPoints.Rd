\name{volTriPoints}
\alias{volTriPoints}
\docType{data}
\title{
   Solute plume concentration data set.
}
\description{
Example solute plume concentration data set
associated with the GWSDAT (\dQuote{GroundWater
Spatiotemporal Data Analysis Tool}) project
\url{https://protect-au.mimecast.com/s/demRC91WzLH6qo3TorzN7?domain=gwsdat.net}.
The \code{deldir} package is used in this project as part of a
numerical routine to estimate plume quantities (mass, average
concentration and centre of mass).
}
\usage{volTriPoints}
\format{
  A data frame with 232 observations on the following 3 variables.
  \describe{
    \item{\code{x}}{The \code{x}-coordinates of the centres of mass
          of the plumes.}
    \item{\code{y}}{The \code{y}-coordinates of the centres of mass
          of the plumes.}
    \item{\code{z}}{The plume concentrations.}
  }
}
\details{
  This data set played a critical role in revealing a bug in the
  Fortran code underlying the \code{deldir()} function.
}
\source{
  These data were kindly provided by Wayne W. Jones of the GWSDAT
  project.  The data set was originally named \code{Vol.Tri.Points};
  the name was changed so as to be more consistent with my usual
  naming conventions.
}
\references{
Jones, W. R., Spence, M. J., Bowman, A. W., Evers, L. and Molinari,
D. A.  (2014). A software tool for the spatiotemporal analysis
and reporting of groundwater monitoring data. \emph{Environmental
Modelling & Software} \bold{55}, pp. 242--249.
}
\examples{
dvtp <- deldir(volTriPoints)
plot(dvtp)
}
\keyword{datasets}
