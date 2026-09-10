\name{wpstCLASS}
\alias{wpstCLASS}
\title{
Predict values using new time series values via a non-decimated wavelet packet
discrimination object.
}
\description{Given a timeseries (\code{timeseries}) and another time series
of categorical values (\code{groups}) the \code{\link{makewpstDO}} produces
a model that permits discrimination of the \code{groups} series using
a discriminant analysis based on a restricted set of non-decimated
wavelet packet coefficients of \code{timeseries}. The current function
enables new \code{timeseries} data, to be used in conjunction with
the model to generate new, predicted, values of the \code{groups} time series.
}
\usage{
wpstCLASS(newTS, wpstDO)
}
\arguments{
\item{newTS}{A new segment of time series values, of the same time series
	that was used as the dependent variable used to construct the
	wpstDO object}
\item{wpstDO}{An object that uses values of a dependent time series to
	build a discriminatory model of a groups time series. Output
	from the \code{\link{makewpstDO}} function}
}
\details{
This function performs the same nondecimated wavelet packet (NDWPT) transform
of the \code{newTS} data that was used to analyse the original \code{timeseries}
and the details of this transform are stored within the \code{wpstDO}
object. Then, using information that was recorded in \code{wpstDO}
the packets with the same level/index are extracted from the new NDWPT 
and formed into a matrix. Then the linear discriminant variables,
again stored in \code{wpstDO} are used to form predictors of the
original \code{groups} time series, ie new values of \code{groups}
that correspond to the new values of \code{timeseries}.
}
\value{
The prediction using the usual R \code{predict.lda} function. The
predicted values are stored in the \code{class} component of that list.
}
\seealso{\code{\link{makewpstDO}}}
\examples{
#
# See example at the end of help page for makewpstDO
#
}
\author{G P Nason}
\keyword{ts}
\keyword{multivariate}
