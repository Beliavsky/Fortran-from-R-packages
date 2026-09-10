\name{threshold.wst}
\alias{threshold.wst}
\title{Threshold (NDWT) packet-ordered non-decimated wavelet decomposition object}
\description{
This function provides various ways to threshold a \code{\link{wst}} class object
}
\usage{
\method{threshold}{wst}(wst, levels = 3:(nlevelsWT(wst) - 1), dev = madmad, policy = 
	"universal", value = 0, by.level = FALSE, type = "soft", verbose
	 = FALSE, return.threshold = FALSE, cvtol = 0.01, cvnorm = l2norm, 
	add.history = TRUE, \dots)
}
\arguments{
\item{wst}{The packet ordered non-decimated wavelet decomposition object that you wish to threshold.}
\item{levels}{a vector of integers which determines which scale levels are thresholded in the decomposition. Each integer in the vector must refer to a valid level in the \code{\link{wst}} object supplied. This is usually any integer from 0 to \code{\link{nlevelsWT}}(wst)-1 inclusive. Only the levels in this vector contribute to the computation of the threshold and its application. }
\item{dev}{this argument supplies the function to be used to compute the spread of the absolute values coefficients. The function supplied must return a value of spread on the variance scale (i.e. not standard deviation) such as the \code{var()} function. A popular, useful and robust alternative is the \code{\link{madmad}} function}
\item{policy}{selects the technique by which the threshold value is selected. Each policy corresponds to a method in the literature. At present the different policies are: "\code{universal}", "\code{LSuniversal}", "\code{\link{sure}}", "\code{cv}", "\code{manual}", The policies are described in detail \code{below}. }
\item{value}{This argument conveys the user supplied threshold. If the \code{policy="manual"} then value is the actual threshold value.}
\item{by.level}{If FALSE then a global threshold is computed on and applied to all scale levels defined in \code{levels}. If TRUE a threshold is computed and applied separately to each scale level.}
\item{type}{determines the type of thresholding this can be "\code{hard}" or "\code{soft}".}
\item{verbose}{if TRUE then the function prints out informative messages as it progresses.} 
\item{return.threshold}{If this option is TRUE then the actual \emph{value} of the threshold is returned. If this option is FALSE then a thresholded version of the input is returned.}
\item{cvtol}{Parameter for the cross-validation "\code{cv}" policy.} 
\item{cvnorm}{A function to compute the distance between two vectors. Two useful possibilities are \code{\link{l2norm}} and \code{\link{linfnorm}}. Selection of different metrics causes the cross-validation denoising method to optimize for different criteria.}
\item{add.history}{If \code{TRUE} then the thresholding operation details are add to the returned \code{\link{wst}} object. This can be useful when later tracing how an object has been treated.} 
\item{\dots}{any other arguments}
}
\details{
This function thresholds or shrinks wavelet coefficients stored in a \code{\link{wst}} object and returns the coefficients in a modified \code{\link{wst}} object. The thresholding step is an essential component of denoising using the \code{packet-ordered non-decimated wavelet transform}. If the denoising is carried out using the \code{\link{AvBasis}} basis averaging technique then this software is an implementation of the Coifman and Donoho translation-invariant (TI) denoising. (Although it is the denoising technique which is translation invariant, not the packet ordered non-decimated transform, which is translation equivariant). However, the \code{threshold.wst} algorithm can be used in other denoising techniques in conjunction with the basis selection and inversion functions \code{\link{MaNoVe}} and \code{\link{InvBasis}}.
 
The basic idea of thresholding is very simple. In a signal plus noise model the wavelet transform of signal is very sparse, the wavelet transform of noise is not (in particular, if the noise is iid Gaussian then so if the noise contained in the wavelet coefficients). Thus since the signal gets concentrated in the wavelet coefficients and the noise remains "spread" out it is "easy" to separate the signal from noise by keeping large coefficients (which correspond to signal) and delete the small ones (which correspond to noise). However, one has to have some idea of the noise level (computed using the dev option in threshold functions). If the noise level is very large then it is possible, as usual, that no signal "sticks up" above the noise. 

Many of the pragmatic comments for successful thresholding given in the help for \code{\link{threshold.wd}} hold true here: after all non-decimated wavelet transforms are merely organized collections of standard (decimated) discrete wavelet transforms. We reproduce some of the issues relevant to thresholding \code{\link{wst}} objects. 

Some issues to watch for: 

\describe{
\item{levels}{The default of \code{levels = 3:(nlevelsWT(wd) - 1)} for the \code{levels} option most certainly does not work globally for all data problems and situations. The level at which thresholding begins (i.e. the given threshold and finer scale wavelets) is called the \code{primary resolution} and is unique to a particular problem. In some ways choice of the primary resolution is very similar to choosing the bandwidth in kernel regression albeit on a logarithmic scale. See Hall and Patil, (1995) and Hall and Nason (1997) for more information. For each data problem you need to work out which is the best primary resolution. This can be done by gaining experience at what works best, or using prior knowledge. It is possible to "automatically" choose a "best" primary resolution using cross-validation (but not yet in WaveThresh). 

Secondly the levels argument computes and applies the threshold at the levels specified in the \code{levels} argument. It does this for all the levels specified. Sometimes, in wavelet shrinkage, the threshold is computed using only the finest scale coefficients (or more precisely the estimate of the overall noise level). If you want your threshold variance estimate only to use the finest scale coefficients (e.g. with universal thresholding) then you will have to apply the \code{threshold.wd} function twice. Once (with levels set equal to \code{\link{nlevelsWT}}(wd)-1 and with \code{return.threshold=TRUE} to return the threshold computed on the finest scale and then apply the threshold function with the \code{manual} option supplying the value of the previously computed threshold as the \code{value} options. }
\item{by.level}{for a \code{\link{wd}} object which has come from data with noise that is correlated then you should have a threshold computed for each resolution level. See the paper by Johnstone and Silverman, 1997.}
}
}
\value{
An object of class \code{\link{wst}}. This object contains the thresholded wavelet coefficients. Note that if the \code{return.threshold} option is set to TRUE then the threshold values will be returned rather than the thresholded object. 
}
\note{
This section gives a brief description of the different thresholding policies available. For further details \emph{see the associated papers}. If there is no paper available then a small description is provided here. More than one policy may be good for problem, so experiment! Some of the policies here were specifically adapted to the This section gives a brief description of the different thresholding policies available. For further details see the associated papers. If there is no paper available then a small description is provided here. More than one policy may be good for problem, so experiment! Some of the policies here were specifically adapted to the \code{\link{wst.object}} but some weren't so beware. They are arranged here in alphabetical order:
 
\describe{
\item{cv}{See Nason, 1996.}
\item{LSuniversal}{See Nason, von Sachs and Kroisandt, 1998. This is used for smoothing of a wavelet periodogram and shouldn't be used generally.}
\item{manual}{specify a user supplied threshold using \code{value} to pass the value of the threshold. The \code{value} argument should be a vector. If it is of length 1 then it is replicated to be the same length as the \code{levels} vector, otherwise it is repeated as many times as is necessary to be the levels vector's length. In this way, different thresholds can be supplied for different levels. Note that the \code{by.level} option has no effect with this policy.}\item{sure}{See Donoho and Johnstone, 1994 and Johnstone and Silverman, 1997.}
\item{universal}{See Donoho and Johnstone, 1995.} 
}
}
\section{RELEASE}{Version 3.6 Copyright Guy Nason 1997 }
\seealso{
\code{\link{AvBasis}}, \code{\link{AvBasis.wst}}, \code{\link{InvBasis}}, \code{\link{InvBasis.wst}}, \code{\link{MaNoVe}},\code{\link{MaNoVe.wst}}, \code{\link{wst}}, \code{\link{wst.object}}, \code{\link{threshold}}. 
}
\keyword{smooth}
\keyword{nonlinear}
\author{G P Nason}
