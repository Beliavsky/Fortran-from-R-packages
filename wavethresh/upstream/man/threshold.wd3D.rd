\name{threshold.wd3D}
\alias{threshold.wd3D}
\title{Threshold 3D DWT object}
\description{
This function provides various ways to threshold a \code{\link{wd3D}} class object. 
}
\usage{
\method{threshold}{wd3D}(wd3D, levels = 3:(nlevelsWT(wd3D) - 1), type = "hard", policy = 
        "universal", by.level = FALSE, value = 0, dev = var, verbose = FALSE, 
        return.threshold = FALSE, \dots)
}
\arguments{
\item{wd3D}{The 3D DWT wavelet decomposition object that you wish to threshold.}
\item{levels}{a vector of integers which determines which scale levels are thresholded in the decomposition. Each integer in the vector must refer to a valid level in the \code{\link{wd3D}} object supplied. This is usually any integer from 0 to \code{\link{nlevelsWT}}(wd3D)-1 inclusive. Only the levels in this vector contribute to the computation of the threshold and its application. }
\item{type}{determines the type of thresholding this can be "\code{hard}" or "\code{soft}".}
\item{policy}{selects the technique by which the threshold value is selected. Each policy corresponds to a method in the literature. At present the different policies are: "\code{universal}" and "\code{manual}". The policies are described in detail \code{below}.} 
\item{by.level}{If FALSE then a global threshold is computed on and applied to all scale levels defined in \code{levels}. If TRUE a threshold is computed and applied separately to each scale level.}
\item{value}{This argument conveys the user supplied threshold. If the \code{policy="manual"} then value is the actual threshold value.} 
\item{dev}{this argument supplies the function to be used to compute the spread of the absolute values coefficients. The function supplied must return a value of spread on the variance scale (i.e. not standard deviation) such as the \code{var()} function. A popular, useful and robust alternative is the \code{\link{madmad}} function.}
\item{verbose}{if TRUE then the function prints out informative messages as it progresses.}
\item{return.threshold}{If this option is TRUE then the actual \emph{value} of the threshold is returned. If this option is FALSE then a thresholded version of the input is returned.} 
\item{\dots}{any other arguments}
}
\details{
This function thresholds or shrinks wavelet coefficients stored in a \code{\link{wd3D}} object and returns the coefficients in a modified \code{\link{wd3D}} object. See the seminal papers by Donoho and Johnstone for explanations about thresholding. For a gentle introduction to wavelet thresholding (or shrinkage as it is sometimes called) see Nason and Silverman, 1994. For more details on each technique see the descriptions of each method below 

The basic idea of thresholding is very simple. In a signal plus noise model the wavelet transform of signal is very sparse, the wavelet transform of noise is not (in particular, if the noise is iid Gaussian then so if the noise contained in the wavelet coefficients). Thus since the signal gets concentrated in the wavelet coefficients and the noise remains "spread" out it is "easy" to separate the signal from noise by keeping large coefficients (which correspond to signal) and delete the small ones (which correspond to noise). However, one has to have some idea of the noise level (computed using the dev option in threshold functions). If the noise level is very large then it is possible, as usual, that no signal "sticks up" above the noise. 

There are many components to a successful thresholding procedure. Some components have a larger effect than others but the effect is not the same in all practical data situations. Here we give some rough practical guidance, although \emph{you must refer to the papers below when using a particular technique.} \bold{You cannot expect to get excellent performance on all signals unless you fully understand the rationale and limitations of each method below.} I am not in favour of the "black-box" approach. The thresholding functions of WaveThresh3 are not a black box: experience and judgement are required! 

Some issues to watch for: 

\describe{
\item{levels}{The default of \code{levels = 3:(wd$nlevelsWT - 1)} for the \code{levels} option most certainly does not work globally for all data problems and situations. The level at which thresholding begins (i.e. the given threshold and finer scale wavelets) is called the \emph{primary resolution} and is unique to a particular problem. In some ways choice of the primary resolution is very similar to choosing the bandwidth in kernel regression albeit on a logarithmic scale. See Hall and Patil, (1995) and Hall and Nason (1997) for more information. For each data problem you need to work out which is the best primary resolution. This can be done by gaining experience at what works best, or using prior knowledge. It is possible to "automatically" choose a "best" primary resolution using cross-validation (but not in WaveThresh). 

Secondly the levels argument computes and applies the threshold at the levels specified in the \code{levels} argument. It does this for all the levels specified. Sometimes, in wavelet shrinkage, the threshold is computed using only the finest scale coefficients (or more precisely the estimate of the overall noise level). If you want your threshold variance estimate only to use the finest scale coefficients (e.g. with universal thresholding) then you will have to apply the \code{threshold.wd} function twice. Once (with levels set equal to \code{\link{nlevelsWT}}(wd)-1 and with \code{return.threshold=TRUE} to return the threshold computed on the finest scale and then apply the threshold function with the manual option supplying the value of the previously computed threshold as the value options. }

\item{by.level}{for a \code{\link{wd}} object which has come from data with noise that is correlated then you should have a threshold computed for each resolution level. See the paper by Johnstone and Silverman, 1997.}
}
}
\value{
An object of class \code{\link{wd3D}}. This object contains the thresholded wavelet coefficients. Note that if the return.threshold option is set to TRUE then the threshold values will be returned rather than the thresholded object. 
}
\note{
POLICIES
 
This section gives a brief description of the different thresholding policies available. For further details \emph{see the associated papers}. If there is no paper available then a small description is provided here. More than one policy may be good for problem, so experiment! They are arranged here in alphabetical order: 

\describe{
\item{manual}{specify a user supplied threshold using value to pass the value of the threshold. The value argument should be a vector. If it is of length 1 then it is replicated to be the same length as the \code{levels} vector, otherwise it is repeated as many times as is necessary to be the \code{levels} vector's length. In this way, different thresholds can be supplied for different levels. Note that the \code{by.level} option has no effect with this policy. }
\item{universal}{See Donoho and Johnstone, 1995.}
}
}
\section{RELEASE}{Version 3.9.6 Copyright Guy Nason 1997.}
\seealso{
\code{\link{threshold}}, \code{\link{accessD.wd3D}}, \code{\link{print.wd3D}}, \code{\link{putD.wd3D}}, \code{\link{putDwd3Dcheck}}, \code{\link{summary.wd3D}}, \code{\link{threshold.wd3D}}, \code{\link{wd3D.object}}, \code{\link{wr3D}}. 
}
\examples{
#
# Generate some test data
#
test.data <- array(rnorm(8*8*8), dim=c(8,8,8))
testwd3D <- wd3D(test.data)
#
# Now let's threshold
#
testwd3DT <- threshold(testwd3D, levels=1:2)
#
# That's it, one can apply wr3D now to reconstruct
# if you like!
#
}
\keyword{smooth}
\keyword{nonlinear}
\author{G P Nason}
