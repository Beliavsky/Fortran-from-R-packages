\name{makewpstDO}
\alias{makewpstDO}
\title{Help page for a function}
\description{Takes two time series: one a real-valued discrete-time time
series, timeseries, the other, groups, a time series containing factor levels. This function
performs a discriminant analysis of groups on a subset of the best-correlating
nondecimated wavelet packets of timeseries
}
\usage{
makewpstDO(timeseries, groups, filter.number=10, family="DaubExPhase",
	mincor=0.69999999999999996)
}
\arguments{
\item{timeseries}{The time series which is the `dependent variable',
	ie discrimination will be performed on the variables extracted
	from the non-decimated wavelet packet transform of this time series}
\item{groups}{The factor levels as a time series}
\item{filter.number}{The smoothness of the wavelet involved in the
	nondecimated wavelet packet transform. See \code{\link{filter.select}}}
\item{family}{The wavelet family, see \code{\link{filter.select}}}
\item{mincor}{Variables from the nondecimated wavelet packet transform
	with correlations less than this argument will be discarded in the
	first pass, and not considered as possible useful discriminants}
}
\details{
This function implements the `discrimination' version of the 
"Wavelet packet transfer function modelling of nonstationary series"
by Guy Nason and Theofanis Sapatinas, \emph{Statistics and Computing},
\bold{12}, 45-56.

The function first takes the non-decimated wavelet packet transform
of \code{timeseries}
using the \code{\link{wpst}} function. Then the set of nondecimated wavelet
packets is put into matrix form using the \code{\link{wpst2discr}}
function. The \code{\link{Best1DCols}} function selects those variables
from the matrix whose correlation with the \code{groups} time series
is greater than \code{mincor}. The selected variables are put into a
reduced matrix. 

The next step, \code{\link{BMdiscr}},
performs a linear discriminant analysis of the
\code{groups} values onto the reduced matrix. In principle, one could have
carried out a discriminant analysis using the full matrix of all the
packets, but the problem is not well-conditioned and computationally
efficient. The strategy adopted by Nason and Sapatinas is to do a
"first pass" to select a large number of "likely" variables that might
contribute something to discrimination, and then carry out a "second pass"
which performs a more detailed analysis to jointly determine which variables
are the key ones for discrimination. 

Note, using the discriminant model developed here, it is possible to
use future values of \code{timeseries} and the model to predict future
values of \code{groups}. See example below.
}
\value{
An object of class \code{wpstDO}. This is a list containing the following
components.

\item{BPd}{Object returned from the \code{\link{BMdiscr}} function.
	Contains the reduced matrix and the discriminant object}
\item{BP}{Object returned from the \code{\link{Best1DCols}} function,
	essentially the reduced matrix and the groups variable.}
\item{filter}{The details of the wavelet filter used. This is used
	if the other components are used to perform discrimination on
	new data one needs to know what wavelet was used to perform
	the original nondecimated wavelet packet transform.} 
}
\seealso{
\code{\link{basisplot.BP}},
\code{\link{Best1DCols}},
\code{\link{BMdiscr}},
\code{\link{wpst}},
\code{\link{wpst2discr}},
\code{\link{wpstCLASS}}
}
\examples{
#
# Use BabySS and BabyECG data for this example.
#
# Want to predict future values of BabySS from future values of BabyECG
#
# Build model on first 256 values of both
#
data(BabyECG)
data(BabySS)
BabyModel <- makewpstDO(timeseries=BabyECG[1:256], groups=BabySS[1:256],
	mincor=0.5)
#
# The results (ie print out answer)
#BabyModel
#Stationary wavelet packet discrimination object
#Composite object containing components:[1] "BPd"    "BP"     "filter"
#Fisher's discrimination: done
#BP component has the following information
#BP class object. Contains "best basis" information
#Components of object:[1] "nlevelsWT"     "BasisMatrix" "level"       "pkt"         "basiscoef"
#[6] "groups"
#Number of levels  8
#List of "best" packets
#Level id Packet id Basis coef
#[1,]        4         0  0.7340580
#[2,]        5         0  0.6811251
#[3,]        6         0  0.6443167
#[4,]        3         0  0.6193434
#[5,]        7         0  0.5967620
#[6,]        0         3  0.5473777
#[7,]        1        53  0.5082849
#
# You can plot the select basis graphically using
#
\dontrun{basisplot(BabyModel$BP)}
#
# An interesting thing are the final "best" packets, these form the
# "reduced" matrix, and the final discrimination is done on this
# In this case 7 wavelet packets were identified as being good for
# univariate high correlation.
#
# In the second pass lda analysis, using the reduced matrix, the following
# turns up as the best linear discriminant vectors
#
# The discriminant variables can be obtained by typing
#BabyModel$BPd$dm$scaling
#LD1        LD2
#[1,] 5.17130434  1.8961807
#[2,] 1.56487144 -3.5025251
#[3,] 1.69328553  1.1585477
#[4,] 3.63362324  8.4543247
#[5,] 0.15202947 -0.4530523
#[6,] 0.35659009 -0.3850318
#[7,] 0.09429836 -0.1281240
#
#
# Now, suppose we get some new data for the BabyECG time series.
# For the purposes of this example, this is just the continuing example
# ie BabyECG[257:512]. We can use our new discriminant model to predict
# new values of BabySS
#
BabySSpred <- wpstCLASS(newTS=BabyECG[257:512], BabyModel)
#
# Let's look at the first 10 (eg) values of this prediction
#
#BabySSpred$class[1:10]
#[1] 4 4 4 4 4 4 4 4 4 4
#Good. Now let's look at what the "truth" was:
#BabySS[257:267]
#[1] 4 4 4 4 4 4 4 4 4 4
#Good. However, the don't agree everywhere, let's do a cross classification
#between the prediction and the truth.
#
#> table(tmp2$class, BabySS[257:512])
#
#      1   2   3   4
#  1   4   1   1   0
#  2 116   0  23   3
#  4   2  12   0  94
#
#So class 3 and 4 agree pretty much, but class 1 has been mispredicted at class
#2 a lot.
}
\author{G P Nason}
\keyword{multivariate}
\keyword{ts}
