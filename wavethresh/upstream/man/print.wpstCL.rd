\name{print.wpstCL}
\alias{print.wpstCL}
\title{Prints some information about a wpstCL object}
\usage{
\method{print}{wpstCL}(x, \dots)
}
\arguments{
\item{x}{wpstCL object to print info about}
\item{\dots}{Other arguments}
}
\description{
Prints basic information about a wpstCL object
}
\details{
Description says all
}
\value{
Nothing
}
\seealso{\code{\link{makewpstDO}},\code{\link{wpstCLASS}}}
\examples{
#
# Use BabySS and BabyECG data for this example.
#
# Want to predict future values of BabySS from future values of BabyECG
#
# Build model on first 256 values of both
#
# See example in makewpstDO from which this one originates
#
data(BabyECG)
data(BabySS)
BabyModel <- makewpstDO(timeseries=BabyECG[1:256], groups=BabySS[1:256],
	mincor=0.5)
#
# Now, suppose we get some new data for the BabyECG time series.
# For the purposes of this example, this is just the continuing example
# ie BabyECG[257:512]. We can use our new discriminant model to predict
# new values of BabySS
#
BabySSpred <- wpstCLASS(newTS=BabyECG[257:512], BabyModel)
#
BabySSpred
#wpstCL class object
#Results of applying discriminator to time series
#Components:  BasisMatrix BasisMatrixDM wpstDO PredictedOP PredictedGroups 
}
\author{G P Nason}
\keyword{print}
