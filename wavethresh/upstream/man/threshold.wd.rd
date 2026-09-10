\name{threshold.wd}
\alias{threshold.wd}
\title{Threshold (DWT) wavelet decomposition object}
\description{
This function provides various ways to threshold a \code{\link{wd}} class object. 
}
\usage{
\method{threshold}{wd}(wd, levels = 3:(nlevelsWT(wd) - 1), type = "soft", policy = "sure",
        by.level = FALSE, value = 0, dev = madmad, boundary = FALSE, verbose = FALSE,
        return.threshold = FALSE, force.sure = FALSE, cvtol = 0.01,
	cvmaxits=500, Q = 0.05, OP1alpha = 0.05, 
        alpha = 0.5, beta = 1, C1 = NA, C2 = NA, C1.start = 100,
	al.check=TRUE, \dots)
}
\arguments{
\item{wd}{The DWT wavelet decomposition object that you wish to threshold.}
\item{levels}{a vector of integers which determines which scale levels are thresholded in the decomposition. Each integer in the vector must refer to a valid level in the \code{\link{wd}} object supplied. This is usually any integer from 0 to \code{\link{nlevelsWT}}(wd)-1 inclusive. Only the levels in this vector contribute to the computation of the threshold and its application.}
\item{type}{determines the type of thresholding this can be "hard" or "soft".}
\item{policy}{selects the technique by which the threshold value is selected. Each policy corresponds to a method in the literature. At present the different policies are: "\code{universal}", "\code{LSuniversal}", "\code{\link{sure}}", "\code{BayesThresh}", "\code{cv}", "\code{fdr}", "\code{op1}", "\code{op2}", "\code{manual}", "\code{mannum}" and "\code{probability}". The policies are described in detail below.}
\item{by.level}{If FALSE then a global threshold is computed on and applied to all scale levels defined in \code{levels}. If TRUE a threshold is computed and applied separately to each scale level.}
\item{value}{This argument conveys the user supplied threshold. If the \code{policy="manual"} then \code{value} is the actual threshold value; if the \code{policy="mannum"} then \code{value} conveys the total number of ordered coefficients kept (from the largest); if \code{policy="probability"} then \code{value} conveys the the user supplied quantile level.}
\item{dev}{this argument supplies the function to be used to compute the spread of the absolute values coefficients. The function supplied must return a value of spread on the variance scale (i.e. not standard deviation) such as the \code{var()} function. A popular, useful and robust alternative is the \code{\link{madmad}} function.}
\item{boundary}{If this argument is TRUE then the boundary bookeeping values are included for thresholding, otherwise they are not.}
\item{verbose}{if TRUE then the function prints out informative messages as it progresses.}
\item{return.threshold}{If this option is TRUE then the actual \emph{value} of the threshold is returned. If this option is FALSE then a thresholded version of the input is returned.}
\item{force.sure}{If TRUE then the \code{\link{sure}} threshold is computed on a vector even when that vector is very sparse. If FALSE then the normal SUREshrink procedure is followed whereby the universal threshold is used for sparse vectors of coefficients.}
\item{cvtol}{Parameter for the cross-validation \code{"cv"} policy.}
\item{cvmaxits}{Maximum number of iterations allowed for the cross-validation \code{"cv"} policy.}
\item{Q}{Parameter for the false discovery rate \code{"fdr"} policy.}
\item{OP1alpha}{Parameter for Ogden and Parzen's first "\code{op1}" and \code{"op2"} policies.}
\item{alpha}{Parameter for BayesThresh \code{"BayesThresh"} policy.}
\item{beta}{Parameter for BayesThresh \code{"BayesThresh"} policy.}
\item{C1}{Parameter for BayesThresh \code{"BayesThresh"} policy.}
\item{C2}{Parameter for BayesThresh \code{"BayesThresh"} policy.} 
\item{C1.start}{Parameter for BayesThresh \code{"BayesThresh"} policy.}
\item{al.check}{If TRUE then the function checks that the levels are
	in ascending order. If they are not then this can be an
	indication that the default level arguments are not appropriate
	for this data set (\code{wd} object). However, a strange order
	might be appropriate for some reason if deliberately set, so setting
	this argument equal to FALSE turns off the check and warning.}
\item{\dots}{any other arguments}
}
\details{
This function thresholds or shrinks wavelet coefficients stored in a \code{\link{wd}} object and returns the coefficients in a modified \code{\link{wd}} object. See the seminal papers by Donoho and Johnstone for explanations about thresholding. For a gentle introduction to wavelet thresholding (or shrinkage as it is sometimes called) see Nason and Silverman, 1994. For more details on each technique see the descriptions of each method below 

The basic idea of thresholding is very simple. In a signal plus noise model the wavelet transform of signal is very sparse, the wavelet transform of noise is not (in particular, if the noise is iid Gaussian then so if the noise contained in the wavelet coefficients). Thus since the signal gets concentrated in the wavelet coefficients and the noise remains "spread" out it is "easy" to separate the signal from noise by keeping large coefficients (which correspond to signal) and delete the small ones (which correspond to noise). However, one has to have some idea of the noise level (computed using the dev option in threshold functions). If the noise level is very large then it is possible, as usual, that no signal "sticks up" above the noise. 

There are many components to a successful thresholding procedure. Some components have a larger effect than others but the effect is not the same in all practical data situations. Here we give some rough practical guidance, although \emph{you must refer to the papers below when using a particular technique.} \bold{You cannot expect to get excellent performance on all signals unless you fully understand the rationale and limitations of each method below.} I am not in favour of the "black-box" approach. The thresholding functions of WaveThresh3 are not a black box: experience and judgement are required! 

Some issues to watch for: 

\describe{
\item{levels}{The default of \code{levels = 3:(wd$nlevelsWT - 1)} for the \code{levels} option most certainly does not work globally for all data problems and situations. The level at which thresholding begins (i.e. the given threshold and finer scale wavelets) is called the \emph{primary resolution} and is unique to a particular problem. In some ways choice of the primary resolution is very similar to choosing the bandwidth in kernel regression albeit on a logarithmic scale. See Hall and Patil, (1995) and Hall and Nason (1997) for more information. For each data problem you need to work out which is the best primary resolution. This can be done by gaining experience at what works best, or using prior knowledge. It is possible to "automatically" choose a "best" primary resolution using cross-validation (but not in WaveThresh). 

Secondly the levels argument computes and applies the threshold at the levels specified in the \code{levels} argument. It does this for all the levels specified. Sometimes, in wavelet shrinkage, the threshold is computed using only the finest scale coefficients (or more precisely the estimate of the overall noise level). If you want your threshold variance estimate only to use the finest scale coefficients (e.g. with universal thresholding) then you will have to apply the \code{threshold.wd} function twice. Once (with levels set equal to \code{\link{nlevelsWT}}(wd)-1 and with \code{return.threshold=TRUE} to return the threshold computed on the finest scale and then apply the threshold function with the manual option supplying the value of the previously computed threshold as the value options.

Thirdly, if you apply wavelet shrinkage to a small data set then you need to ensure you've chosen the \code{levels} argument appropriately. For example,
if your original data was of length 8, then the associated \code{wd}
wavelet decomposition object will only have levels 0, 1 and 2. So, the default
argument for levels (starting at 3 and higher) will almost certainly
be wrong. The code now warns for these situations.
} 


\item{by.level}{for a \code{\link{wd}} object which has come from data with noise that is correlated then you should have a threshold computed for each resolution level. See the paper by Johnstone and Silverman, 1997.} 
}
}
\value{
An object of class \code{\link{wd}}. This object contains the thresholded wavelet coefficients. Note that if the \code{return.threshold} option is set to TRUE then the threshold values will be returned rather than the thresholded object. 
}
\references{Various code segments detailed above were kindly donated by Felix Abramovich, Theofanis Sapatinas and Todd Ogden. }
\note{

POLICIES 
This section gives a brief description of the different thresholding policies available. For further details see \emph{the associated papers}. If there is no paper available then a small description is provided here. More than one policy may be good for problem, so experiment! They are arranged here in alphabetical order: 
\describe{
\item{BayesThresh}{See Abramovich, Silverman and Sapatinas, (1998). Contributed by Felix Abramovich and Fanis Sapatinas.}
\item{cv}{See Nason, 1996.}
\item{fdr}{See Abramovich and Benjamini, 1996. Contributed by Felix Abramovich.}
\item{LSuniversal}{See Nason, von Sachs and Kroisandt, 1998. This is used for smoothing of a wavelet periodogram and shouldn't be used generally.}
\item{manual}{specify a user supplied threshold using \code{value} to pass the value of the threshold. The \code{value} argument should be a vector. If it is of length 1 then it is replicated to be the same length as the \code{levels} vector, otherwise it is repeated as many times as is necessary to be the \code{levels} vector's length. In this way, different thresholds can be supplied for different levels. Note that the \code{by.level} option has no effect with this policy.}
\item{mannum}{You decided how many of the largest (in absolute value) coefficients that you want to keep and supply this number in value.}
\item{op1}{See Ogden and Parzen, 1996. Contributed by Todd Ogden.}
\item{op2}{See Ogden and Parzen, 1996. Contributed by Todd Ogden.}
\item{probability}{The \code{probability} policy works as follows. All coefficients that are smaller than the valueth quantile of the coefficients are set to zero. If \code{by.level} is false, then the quantile is computed for all coefficients in the levels specified by the "levels" vector; if \code{by.level} is true, then each level's quantile is estimated separately. The probability policy is pretty stupid - do not use it.}
\item{sure}{See Donoho and Johnstone, 1994.}
\item{universal}{See Donoho and Johnstone, 1995.} 
}
}
\section{RELEASE}{Version 3.6 Copyright Guy Nason and others 1997 }
\seealso{
\code{\link{wd}}, \code{\link{wd.object}}, \code{\link{wr}}, \code{\link{wr.wd}}, \code{\link{threshold}}. 
}
\examples{
#
# Generate some test data
#
test.data <- example.1()$y
\dontrun{ts.plot(test.data)}
#
# Generate some noisy data
#
ynoise <- test.data + rnorm(512, sd=0.1)
#
# Plot it
#
\dontrun{ts.plot(ynoise)}
#
# Now take the discrete wavelet transform
# N.b. I have no idea if the default wavelets here are appropriate for
# this particular examples. 
#
ynwd <- wd(ynoise)
\dontrun{plot(ynwd)}
#
# Now do thresholding. We'll use a universal policy, 
# and madmad deviance estimate on the finest
# coefficients and return the threshold. We'll also get it to be verbose
# so we can watch the process.
#
ynwdT1 <- threshold(ynwd, policy="universal", dev=madmad,
		levels= nlevelsWT(ynwd)-1, return.threshold=TRUE,
	verbose=TRUE)
# threshold.wd:
# Argument checking
# Universal policy...All levels at once
# Global threshold is:  0.328410967430135 
#
# Why is this the threshold? Well in this case n=512 so sqrt(2*log(n)),
# the universal threshold,
# is equal to 3.53223. Since the noise is about 0.1 (because that's what
# we generated it to be) the threshold is about 0.353.
#
# Now let's apply this threshold to all levels in the noisy wavelet object
#
ynwdT1obj <- threshold(ynwd, policy="manual", value=ynwdT1,
	levels=0:(nlevelsWT(ynwd)-1))
#
# And let's plot it
#
\dontrun{plot(ynwdT1obj)}
#
# You'll see that a lot of coefficients have been set to zero, or shrunk.
#
# Let's try a Bayesian examples this time!
#
ynwdT2obj <- threshold(ynwd, policy="BayesThresh")
#
# And plot the coefficients
#
\dontrun{plot(ynwdT2obj)}
#
# Let us now see what the actual estimates look like
#
ywr1 <- wr(ynwdT1obj)
ywr2 <- wr(ynwdT2obj)
#
# Here's the estimate using universal thresholding
#
\dontrun{ts.plot(ywr1)}
#
# Here's the estimate using BayesThresh 
#
\dontrun{ts.plot(ywr2)}
}
\keyword{smooth}
\keyword{nonlinear}
\author{G P Nason}
