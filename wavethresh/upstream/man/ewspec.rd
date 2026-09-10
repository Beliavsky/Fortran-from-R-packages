\name{ewspec}
\alias{ewspec}
\title{Compute evolutionary wavelet spectrum estimate.}
\description{
This function computes the evolutionary wavelet spectrum (EWS) estimate from a time series (or non-decimated wavelet transform of a time series). The estimate is computed by taking the non-decimated wavelet transform of the time series data, taking its modulus; smoothing using TI-wavelet shrinkage and then correction for the redundancy caused by use of the non-decimated wavelet transform. Options below beginning with smooth. are passed directly to the TI-wavelet shrinkage routines. 
}
\usage{
ewspec(x, filter.number = 10, family = "DaubLeAsymm",
    UseLocalSpec = TRUE, DoSWT = TRUE, WPsmooth = TRUE,
    verbose = FALSE, smooth.filter.number = 10,
    smooth.family = "DaubLeAsymm",
    smooth.levels = 3:(nlevelsWT(WPwst) - 1), smooth.dev = madmad,
    smooth.policy = "LSuniversal", smooth.value = 0, smooth.by.level = FALSE,
    smooth.type = "soft", smooth.verbose = FALSE,
    smooth.cvtol = 0.01, smooth.cvnorm = l2norm,
    smooth.transform = I, smooth.inverse = I)
}
\arguments{
\item{x}{The time series that you want to analyze. (See DETAILS below on how to supply preprocessed versions of the time series which bypass early parts of the ewspec function). }
\item{filter.number }{This selects the index of the wavelet used in the analysis of the time series (i.e. the wavelet basis functions used to model the time series). For Daubechies compactly supported wavelets the filter number is the number of vanishing moments. }
\item{family }{This selects the wavelet family to use in the analysis of the time series (i.e. which wavelet family to use to model the time series). Only use the Daubechies compactly supported wavelets \code{DaubExPhase} and \code{DaubLeAsymm}. }
\item{UseLocalSpec }{If you input a time series for \code{x} then this argument should always be \code{T}. (However, you can precompute the modulus of the non-decimated wavelet transform yourself and supply it as \code{x} in which case the \code{\link{LocalSpec}} call within this function is not necessary and you can set UseLocalSpec equal to \code{F}). }
\item{DoSWT }{If you input a time series for \code{x} then this argument should always be \code{T}. (However, you can precompute the non-decimated wavelet transform yourself and supply it as \code{x} in which case the \code{wd} call within the function will not be necessary and you can set DoSWT equal to \code{F}). }
\item{WPsmooth}{Normally a wavelet periodogram is smoothed before it is corrected. Use \code{WPsmooth=F} is you do not want any wavelet periodogram smoothing (correction is still done). }
\item{verbose}{If this option is \code{T} then informative messages are printed as the function progresses. }
\item{smooth.filter.number }{This selects the index number of the wavelet that smooths each scale of the wavelet periodogram. See \code{\link{filter.select}} for further details on which wavelets you can use. Generally speaking it is a good idea to use a smoother wavelet for smoothing than the one you used for analysis (above) but since one still wants local smoothing it is best not to use a wavelet that is much smoother. }
\item{smooth.family}{This selects the wavelet family that smooths each scale of the wavelet periodogram. See \code{filter.select} for further details on which wavelets you can use. There is no need to use the same family as you used to analyse the time series. }
\item{smooth.levels }{The levels to smooth when performing the TI-wavelet shrinkage smoothing. }
\item{smooth.dev }{The method for estimating the variance of the empirical wavelet coefficients for smoothing purposes. }
\item{smooth.policy }{The recipe for smoothing: determines how the threshold is chosen. See \code{\link{threshold}} for TI-smoothing and choice of potential policies. For EWS estimation \code{LSuniversal} is recommended for thi Chi-squared nature of the periodogram coefficients. However, if the coefficients are transformed (using \code{smooth.transform} and \code{smooth.inverse}) then other, more standard, policies may be appropriate. }
\item{smooth.value }{When a manual policy is being used this argument is used to supply a threshold value. See \code{threshold} for more information. }
\item{smooth.by.level }{If \code{TRUE} then the wavelet shrinkage is performed by computing and applying a separate threshold to each level in the non-decimated wavelet transform of each scale. Note that each scale in the EWS is smoothed separately and independently: and each smooth consists of taking the (second-stage) non-decimated wavelet transform and applying a threshold to each level of a wavelet transformed scale. 

If \code{FALSE} then the same threshold is applied to the non-decimated wavelet transform of a scale. Different thresholds may be computed for different scales (in the time series model) but the threshold will be the same for each level arising from the non-decimated transform of a scale. 

Note: a \code{scale} refers to a set of coefficients coming from a particular scale of the non-decimated wavelet transform of the time series data that \code{models} the time series. A \code{level} refers to the levels of wavelet coefficients obtained from taking the non-decimated wavelet transform of a particular scale.}
\item{smooth.type }{The type of shrinkage: either "hard" or "soft". }
\item{smooth.verbose }{If \code{T} then informative messages concerning the TI-transform wavelet shrinkage are printed.}
\item{smooth.cvtol }{If cross-validated wavelet shrinkage (\code{smooth.policy="cv"}) is used then this argument supplies the cross-validation tolerance. }
\item{smooth.cvnorm}{no description for object}
\item{smooth.transform }{The transform function to use to transform the wavelet periodogram estimate. The wavelet periodogram coefficients are typically chi-squared in nature, a \code{log} transform can pull the coefficients towards normality so that a \code{smooth.policy} for Gaussian data could be used (e.g. \code{universal}). }
\item{smooth.inverse}{the inverse transform of \code{smooth.transform}. }
}
\details{
This function computes an estimate of the evolutionary wavelet spectrum of a time series according to the paper by Nason, von Sachs and Kroisandt. The function works as follows: 
\describe{
\item{1}{The non-decimated wavelet transform of the series is computed.} 
\item{2}{The squared modulus of the non-decimated wavelet transform is computed (this is the raw wavelet periodogram, which is returned).} 
\item{3}{The squared modulus is smoothed using TI-wavelet shrinkage.} 
\item{4}{The smoothed coefficients are corrected using the inverse of the inner product matrix of the discrete non-decimated autocorrelation wavelets (produced using the ipndacw function).} 
}

To display the EWS use the \code{plot}function on the \code{S} component, see the examples below. 

It is possible to supply the non-decimated wavelet transform of the time series and set \code{DoSWT=F} or to supply the squared modulus of the non-decimated wavelet transform using \code{\link{LocalSpec}} and setting \code{UseLocalSpec=F}. This facility saves time because the function is then only used for smoothing and correction. 

}
\value{
A list with the following components: 

\item{S}{The evolutionary wavelet spectral estimate of the input \code{x}. This object is of class \code{\link{wd}} and so can be plotted, printed in the usual way. }
\item{WavPer}{The raw wavelet periodogram of the input \code{x}. The EWS estimate (above) is the smoothed corrected version of the wavelet periodgram. The wavelet periodogram is of class \code{\link{wd}} and so can be plotted, printed in the usual way. }
\item{rm}{This is the matrix A from the paper by Nason, von Sachs and Kroisandt. Its inverse is used to correct the raw wavelet periodogram. This matrix is computed using the \code{\link{ipndacw}} function. }
\item{irm}{The inverse of the matrix A from the paper by Nason, von Sachs and Kroisandt. It is used to correct the raw wavelet periodogram.} 
}
\references{Nason, G.P., von Sachs, R. and Kroisandt, G. (1998). Wavelet processes and adaptive estimation of the evolutionary wavelet spectrum. \emph{Technical Report}, Department of Mathematics University of Bristol/ Fachbereich Mathematik, Kaiserslautern. }
\section{RELEASE}{Version 3.9 Copyright Guy Nason 1998 }
\seealso{
\code{Baby Data}, \code{\link{filter.select}}, \code{\link{ipndacw}}, \code{\link{LocalSpec}}, \code{\link{threshold}} \code{\link{wd}} \code{\link{wd.object}} 
}
\examples{
#
# Apply the EWS estimate function to the baby data
#
}
\keyword{manip}
\author{G P Nason}
