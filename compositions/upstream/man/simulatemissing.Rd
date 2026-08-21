\name{simulatemissings}
\alias{simulateMissings}
\alias{observeWithAdditiveError}
\alias{observeWithDetectionLimit}
\alias{observeWithDetectionlimit}
\title{Artifical simulation of various kinds of missings/polluted data}
\description{
  These are simulation mechanisms to check that missing techniques perform 
  in sensible ways. They just generate additional missings of the various types 
  in a given dataset, according to a specific process.
}
\usage{
simulateMissings(x, dl=NULL, knownlimit=FALSE,
     MARprob=0.0, MNARprob=0.0, mnarity=0.5, SZprob=0.0)
observeWithAdditiveError(x, sigma=dl/dlf, dl=sigma*dlf, dlf=3,
     keepObs=FALSE, digits=NA, obsScale=1,
     class="acomp")
}
\arguments{
  \item{x}{a dataset that should get the missings}
  \item{dl}{the detection limit described in
    \code{\link{clo}}, to impose an artificial detection limit}
  \item{knownlimit}{a boolean indicating wether the actual detection
    limit is still known in the dataset.}
  \item{MARprob}{the probability of occurence of 'Missings At Random' values}
  \item{MNARprob}{the probability of occurrence of 'Missings Not At Random'. 
    The tendency is that small values have a higher probability to
    be missed.}
  \item{mnarity}{a number between 0 and 1 giving the strength of the
    influence of the actual value in becoming a MNAR. 0 means a MAR
    like behavior and 1 means that it is just the smallest values that
    is lost}
  \item{SZprob}{the probability to obtain a structural zero. This is done
    at random like a MAR.}
  \item{sigma}{the standard deviation of the normal distributed extra additive error}
  \item{dlf}{the distance from 0 at which a datum will be considered BDL}
  \item{keepObs}{should the (closed) data without additive error be returned as an attribute?}
  \item{digits}{rounding to be applied to the data with additive error (see Details)}
  \item{obsScale}{rounding to be applied to the data with additive error 
  (see Details). Should be a power of 10.}
  \item{class}{class of the output object}
}
\value{
  A dataset like \code{x} but with some additional missings.

}
\details{
  Without any additional parameters no missings are generated. 
  The procedure to generate MNAR affects all variables.
  
  Function "simulateMissings" is a multipurpose simulator, where each class 
  of missing value is treated separately, and where detection limits are 
  specified as thresholds.
  
  Function "observeWithAdditiveError" simulates data within a very specific 
  framework, where an additive error of sd=\code{sigma} is added to the input data
  \code{x}, and BDLs are generated if a datum is less than \code{dfl} times
  \code{sigma}. Afterwards, the resulting data are rounded as
  \code{round(data/obsScale,digits)*obsScale}, i.e. a certain observation scale
  \code{obsScale} is chosen, and at that scale, only some \code{digits} are kept.
  This framework is typical of chemical analyses, and it generates both BDLs and
  pollution/rounding of (apparently) "right" data.
}
\references{
  van den Boogaart, K., R. Tolosana-Delgado, and M. Bren (2011). The Compositional
  Meaning of a Detection Limit. In Proceedings of the 4th International Workshop on
  Compositional Data Analysis (2011).

  van den Boogaart, K.G., R. Tolosana-Delgado and M. Templ (2014) Regression with
  compositional response having unobserved components or below detection limit 
  values. Statistical Modelling (in press).

  See \link{compositions.missings} for more details.
}
\seealso{
  \link{compositions.missings}
}
\examples{
data(SimulatedAmounts)
x <- acomp(sa.lognormals)
xnew <- simulateMissings(x,dl=0.05,MAR=0.05,MNAR=0.05,SZ=0.05)
acomp(xnew)
plot(missingSummary(xnew))
}
\author{K.Gerald van den Boogaart}
\keyword{datagen}
