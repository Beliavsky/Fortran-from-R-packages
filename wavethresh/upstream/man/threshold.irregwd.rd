\name{threshold.irregwd}
\alias{threshold.irregwd}
\title{hold irregularly spaced wavelet decomposition object}
\description{
This function provides various ways to threshold a \code{\link{irregwd}} class object. 
}
\usage{\method{threshold}{irregwd}(irregwd,
    levels = 3:(nlevelsWT(wd) - 1), type = "hard", policy = "universal",
    by.level = FALSE, value = 0, dev = var, boundary = FALSE,
    verbose = FALSE, return.threshold = FALSE, 
    force.sure=FALSE, cvtol = 0.01, Q = 0.05, alpha=0.05, \dots)
}
\arguments{
\item{irregwd}{The irregularly spaced wavelet decomposition object that you wish to threshold.}
\item{levels}{a vector of integers which determines which scale levels are thresholded in the decomposition. Each integer in the vector must refer to a valid level in the \code{\link{irregwd}} object supplied. This is usually any integer from 0 to \code{\link{nlevelsWT}}(irregwd)-1 inclusive. Only the levels in this vector contribute to the computation of the threshold and its application.}
\item{type}{determines the type of thresholding this can be "hard" or "soft".}
\item{policy}{selects the technique by which the threshold value is selected. Each policy corresponds to a method in the literature. At present the different policies are: \code{"universal"}, \code{"LSuniversal"}, \code{"sure"}, \code{"cv"}, \code{"fdr"}, \code{"op1"}, \code{"op2"}, \code{"manual"}, \code{"mannum"}, \code{"probability"}. A description of the policies can be obtained by clicking on the above links.}
\item{by.level}{If \code{FALSE} then a global threshold is computed on and applied to all scale levels defined in levels. If \code{TRUE} a threshold is computed and applied separately to each scale level.}
\item{value}{This argument conveys the user supplied threshold. If the \code{policy="manual"} then value is the actual threshold value.}
\item{dev}{this argument supplies the function to be used to compute the spread of the absolute values coefficients. The function supplied must return a value of spread on the variance scale (i.e. not standard deviation) such as the \code{var()} function. A popular, useful and robust alternative is the \code{\link{madmad}} function.}
\item{boundary}{If this argument is \code{TRUE} then the boundary bookeeping values are included for thresholding, otherwise they are not.}
\item{verbose}{if \code{TRUE} then the function prints out informative messages as it progresses.}
\item{return.threshold}{If this option is \code{TRUE} then the actual \emph{value} of the threshold is returned. If this option is \code{FALSE} then a thresholded version of the input is returned.}
\item{force.sure}{If \code{TRUE} then the \code{SURE} threshold is computed on a vector even when that vector is very sparse. If \code{FALSE} then the normal SUREshrink procedure is followed whereby the universal threshold is used for sparse vectors of coefficients.}
\item{cvtol}{Parameter for the cross-validation \code{"cv"} policy.}
\item{Q}{Parameter for the false discovery rate \code{"fdr"} policy.}
\item{alpha}{Parameter for Ogden and Parzen's first \code{"op1"} and \code{"op2"} policies.}
\item{\dots}{other arguments}
}
\details{
This function thresholds or shrinks wavelet coefficients stored in a \code{\link{irregwd}} object and returns the coefficients in a modified \code{\link{irregwd}} object. The thresholding step is an essential component of denoising. 

The basic idea of thresholding is very simple. In a signal plus noise model the wavelet transform of signal is very sparse, the wavelet transform of noise is not (in particular, if the noise is iid Gaussian then so if the noise contained in the wavelet coefficients). Thus since the signal gets concentrated in the wavelet coefficients and the noise remains "spread" out it is "easy" to separate the signal from noise by keeping large coefficients (which correspond to signal) and delete the small ones (which correspond to noise). However, one has to have some idea of the noise level (computed using the dev option in threshold functions). If the noise level is very large then it is possible, as usual, that no signal "sticks up" above the noise. 

For thresholding of an \emph{irregularly spaced wavelet decomposition} things are a little different. The original data are irregularly spaced (i.e. [x,y] where the \eqn{x_i} are irregularly spaced) and even if one assumes iid error on the original data once this has been interpolated to a grid by the \code{\link{makegrid}} function the interpolated data values are not independent. The \code{\link{irregwd}} function computes the wavelet transform of the interpolated data but also computes the variance of each coefficient using a fast transform. This variance information is stored in the c component of \code{\link{irregwd}} objects and this function, \code{threshold.irregwd}, makes use of this variance information when thresholding each coefficient. For more details see Kovac and Silverman, 2000 

Some issues to watch for: 

\describe{
\item{levels}{The default of \code{levels = 3:(wd$nlevelsWT - 1)} for the \code{levels} option most certainly does not work globally for all data problems and situations. The level at which thresholding begins (i.e. the given threshold and finer scale wavelets) is called the \emph{primary resolution} and is unique to a particular problem.
In some ways choice of the primary resolution is very similar to choosing the bandwidth in kernel regression albeit on a logarithmic scale. See Hall and Patil, (1995) and Hall and Nason (1997) for more information. For each data problem you need to work out which is the best primary resolution. This can be done by gaining experience at what works best, or using prior knowledge. It is possible to "automatically" choose a "best" primary resolution using cross-validation (but not yet in WaveThresh).
 
Secondly the levels argument computes and applies the threshold at the levels specified in the levels argument. It does this for all the levels specified. Sometimes, in wavelet shrinkage, the threshold is computed using only the finest scale coefficients (or more precisely the estimate of the overall noise level).
If you want your threshold variance estimate only to use the finest scale coefficients (e.g. with universal thresholding) then you will have to apply the \code{threshold.wd} function twice. Once (with levels set equal to \code{nlevelsWT}(wd)-1 and with \code{return.threshold=TRUE} to return the threshold computed on the finest scale and then apply the threshold function with the manual option supplying the value of the previously computed threshold as the value options.}
\item{by.level}{for a \code{\link{wd}} object which has come from data with noise that is correlated then you should have a threshold computed for each resolution level. See the paper by Johnstone and Silverman, 1997.}
}

}
\value{
An object of class \code{\link{irregwd}}. This object contains the thresholded wavelet coefficients. Note that if the \code{return.threshold} option is set to \code{TRUE} then the threshold values will be returned rather than the thresholded object. 
}
\section{RELEASE}{Version 3.6 Copyright Guy Nason 1997 }
\seealso{
\code{\link{makegrid}}, \code{\link{irregwd}}, \code{\link{irregwd}} object, \code{\link{accessc}},
}
\examples{
#
# See main examples of these functions in the help to makegrid
#
} 
\keyword{manip}
\author{Arne Kovac}
