\name{LocalSpec.wd}
\alias{LocalSpec.wd}
\title{Compute Nason and Silverman raw or smoothed wavelet periodogram. }
\description{
\emph{This smoothing in this function is now obsolete}. You should now use the function \code{\link{ewspec}}.
 
This function computes the Nason and Silverman raw or smoothed wavelet periodogram as described by Nason and Silverman (1995). 
}
\usage{
\method{LocalSpec}{wd}(wdS, lsmooth="none", nlsmooth=FALSE, prefilter=TRUE,
  verbose=FALSE, lw.number=wdS$filter$filter.number,
  lw.family=wdS$filter$family, nlw.number=wdS$filter$filter.number,
  nlw.family=wdS$filter$family, nlw.policy="LSuniversal",
  nlw.levels=0:(nlevelsWT(wdS) - 1), nlw.type="hard", nlw.by.level=FALSE,
  nlw.value=0, nlw.dev=var, nlw.boundary=FALSE, nlw.verbose=FALSE,
  nlw.cvtol=0.01, nlw.Q=0.05, nlw.alpha=0.05, nlw.transform=I,
  nlw.inverse=I, debug.spectrum=FALSE, \dots)
}
\arguments{
Note that all options beginning "nlw" are only used if nlsmooth=T, i.e. iff NONLINEAR wavelet smoothing is used. 

\item{wdS}{The stationary wavelet transform object that you want to smooth or square.} 
\item{lsmooth}{Controls the LINEAR smoothing. There are three options: "none", "Fourier" and "wavelet". They are described below. Note that Fourier begins with a capital "F".}
\item{nlsmooth}{A switch to turn on (or off) the NONLINEAR wavelet shrinkage of (possibly LINEAR smoothed) local power coefficients. This option is either TRUE (to turn on the smoothing) or FALSE (to turn it off).} 
\item{prefilter}{If TRUE then apply a prefilter to the actual stationary wavelet coefficients at each level. This is a low-pass filter that cuts off all frequencies above the highest frequency allowed by the (Littlewood-Paley) wavelet that bandpassed the current level coefficients. If FALSE then no prefilter is applied.}
\item{verbose}{If TRUE then the function chats about what it is doing. Otherwise it is silent.} 
\item{lw.number}{If wavelet LINEAR smoothing is used then this option controls the \code{filter number} of the wavelet within the family used to perform the LINEAR wavelet smoothing.} 
\item{lw.family}{If wavelet LINEAR smoothing is used then this option controls the \code{\link{family}} of the wavelet used to perform the LINEAR wavelet smoothing.}
\item{nlw.number}{If NONLINEAR wavelet smoothing is also used then this option controls the \code{filter number} of the wavelet used to perform the wavelet shrinkage.} 
\item{nlw.family}{If NONLINEAR wavelet smoothing is also used then this option controls the \code{\link{family}} of the wavelet used to perform the wavelet shrinkage.} 
\item{nlw.policy}{If NONLINEAR wavelet smoothing is also used then this option controls the levels to use when performing wavelet shrinkage (see \code{\link{threshold.wd}} for different policy choices).} 
\item{nlw.levels}{If NONLINEAR wavelet smoothing is also used then this option controls the levels to use when performing wavelet shrinkage (see \code{\link{threshold.wd}} for a detailed description of how levels can be chosen).} 
\item{nlw.type}{If NONLINEAR wavelet smoothing is also used then this option controls the type of thresholding used in the wavelet shrinkage (either "hard" or "soft", but see \code{\link{threshold.wd}} for a list). }
\item{nlw.by.level}{If NONLINEAR wavelet smoothing is also used then this option controls whether level-by-level thresholding is used or if one threshold is chosen for all levels (see \code{\link{threshold.wd}}).} 
\item{nlw.value}{If NONLINEAR wavelet smoothing is also used then this option controls if a manual (or similar) policy is supplied to \code{nlw.policy} then the nlw.value option carries the manual threshold value (see \code{\link{threshold.wd}}). }
\item{nlw.dev}{If NONLINEAR wavelet smoothing is also used then this option controls the type of variance estimator that is used in wavelet shrinkages (see \code{\link{threshold.wd}}). One possibility is the Splus var() function, another is the WaveThresh function \code{\link{madmad}}().} 
\item{nlw.boundary}{If NONLINEAR wavelet smoothing is also used then this option controls whether boundary coefficients are also thresholded (see \code{\link{threshold.wd}}).} 
\item{nlw.verbose}{If NONLINEAR wavelet smoothing is also used then this option controls whether the threshold function prints out messages as it thresholds levels (see \code{\link{threshold.wd}}).} 
\item{nlw.cvtol}{If NONLINEAR wavelet smoothing is also used then this option controls the optimization tolerance is cross-validation wavelet shrinkage is used (see \code{\link{threshold.wd}})}
\item{nlw.Q}{If NONLINEAR wavelet smoothing is also used then this option controls the Q value for wavelet shrinkage (see \code{\link{threshold.wd}}).} 
\item{nlw.alpha}{If NONLINEAR wavelet smoothing is also used then this option controls the alpha value for wavelet shrinkage (see \code{\link{threshold.wd}}).} 
\item{nlw.transform}{If NONLINEAR wavelet smoothing is also used then this option controls a transformation that is applied to the squared (and possibly linear smoothed) stationary wavelet coefficients before shrinkage. So, for examples, you might want to set \code{nlw.transform=log} to perform wavelet shrinkage on the logs of the squared (and possibly linear smoothed) stationary wavelet coefficients. }
\item{nlw.inverse}{If NONLINEAR wavelet smoothing is also used then this option controls the inverse transformation that is applied to the wavelet shrunk coefficients before they are put back into the stationary wavelet transform structure. So, for examples, if the \code{nlw.transform} is \code{log()} you should set the inverse to \code{nlw.inverse=exp}.} 
\item{debug.spectrum}{If this option is \code{T} then spectrum plots are produced at each stage of the squaring/smoothing. Therefore if you put in the non-decimated wavelet transform of white noise you can get a fair idea of how the coefficients are filtered at each stage.}
\item{\dots}{any other arguments}
}
\details{
\emph{This smoothing in this function is now obsolete.} Use the function \code{\link{ewspec}} instead. However, this function is still useful for computing the raw periodogram. 

This function attempts to produce a picture of local time-scale power of a signal. There are two main components to this function: linear smoothing of squared coefficients and non-linear smoothing of these. Neither, either or both of these components may be used to process the data. The function expects a non-decimated wavelet transform object (of class wd, type="station") such as that produced by the \code{\link{wd}}() function with the type option set to "\code{station}". The following paragraphs describe the various methods of smoothing. 

\bold{LINEAR SMOOTHING}. There are three varieties of linear smoothing. None simply squares the coefficients. Fourier and wavelet apply linear smoothing methods in accordance to the prescription given in Nason and Silverman (1995). Each level in the SWT corresponds to a band-pass filtering to a frequency range [sl, sh]. After squaring we obtain power in the range [0, 2sl] and [2sl, 2sh]. The linear smoothing gets rid of the power in [2sl, 2sh]. The Fourier method simply applies a discrete Fourier transform (rfft) and cuts off frequencies above 2sl. The wavelet method is a bit more suble. The DISCRETE wavelet transform is taken of a level (i) and all levels within the DWT, j, where j>i are set to zero and then the inverse is taken. Approximately this performs the same operation as the Fourier method only faster. By default the same wavelets are used to perform the linear smoothing as were used to compute the stationary wavelet transform in the first place. This can be changed by altering \code{lw.number} and \code{lw.family}.  

\bold{NONLINEAR SMOOTHING}. After either of the linear smoothing options above it is possible to use wavelet shrinkage upon each level in the squared (and possibly Fourier or wavelet linear smoothed) to denoise the coefficients. This process is akin to smoothing the ordinary periodogram. All the usual wavelet shrinkage options are available as \code{nlw}.* where * is one of the usual \code{\link{threshold.wd}} options. By default the same wavelets are used to perform the wavelet shrinkage as were used to compute the non-decimated wavelet transform. These wavelets can be replaced by altering \code{nlw.number} and \code{nlw.family}. Also, it is possible to transform the squared (and possibly smoothed coefficients) before applying wavelet shrinkage. The transformation is effected by supplying an appropriate transformation function (AND ITS INVERSE) to \code{nlw.transform} and \code{nlw.inverse}. (For examples, \code{nlw.transform=log} and\code{ nlw.inverse=exp} might be a good idea). 
}
\value{
An object of class \code{\link{wd}} a time-ordered non-decimated wavelet transform. Each level of the returned object contains a smoothed wavelet periodogram. Note that this is \bold{not} the \emph{corrected} smoothed wavelet periodogram, or the \emph{evolutionary wavelet spectrum}. Use the function \code{\link{ewspec}} to compute the evolutionary wavelet spectrum. 
}
\references{Nason and Silverman, (1995). }
\section{RELEASE}{Version 3.9 Copyright Guy Nason 1998 }
\seealso{
\code{\link{ewspec}}, 
}
\examples{
#
# This function is obsolete. See ewspec()
#
# Compute the raw periodogram of the BabyECG
# data using the Daubechies least-asymmetric wavelet $N=10$.
#
data(BabyECG)
babywdS <- wd(BabyECG, filter.number=10, family="DaubLeAsymm", type="station")
babyWP <- LocalSpec(babywdS, lsmooth = "none", nlsmooth = FALSE)
\dontrun{plot(babyWP, main="Raw Wavelet Periodogram of Baby ECG")}
#
# Note that the lower levels of this plot are too large. This is partly because
# there are "too many" coefficients at the lower levels. For a better
# picture of the local spectral properties of this time series see
# the examples section of ewspec
#
# Other results of this function can be seen in the paper by
# Nason and Silverman (1995) above.
#
}
\keyword{smooth}
\author{G P Nason}
