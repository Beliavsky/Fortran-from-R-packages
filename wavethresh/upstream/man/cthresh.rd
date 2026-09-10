\name{cthresh}
\alias{cthresh}
\title{Estimate real signal using complex-valued wavelets }
\description{
Implements the multiwavelet style and empirical Bayes shrinkage procedures described in Barber & Nason (2004) 
}
\usage{
cthresh(data, j0 = 3, dwwt, dev = madmad, rule = "hard", 
    filter.number = 3.1, family = "LinaMayrand", plotfn = FALSE,  
    TI = FALSE, details = FALSE, policy = "mws", code = "NAG", tol = 0.01)
}
\arguments{
\item{data}{The data to be analysed. This should be real-valued and of length a power of two.}
\item{j0}{Primary resolution level; no thresholding is done below this level.}
\item{dwwt}{description to come}
\item{dev}{A function to be used to estimate the noise level of the data. The function supplied must return a value of spread on the variance scale (i.e. not standard deviation) such as the var() function. A popular, useful and robust alternative is the madmad function.}
\item{rule}{The type of thresholding done. If policy = "mws", available rules are "hard" or "soft"; if policy = "ebayes", then rule can be "hard", "soft" or "mean".}
\item{filter.number, family}{These parameters specify the wavelet used. See \code{\link{filter.select}} for details. 

Also, if filter.number = 5, estimation is done with all the complex-valued wavelets with 5 vanishing moments and the results averaged. If filter.number = 0, then he averaging is over all available complex-valued wavelets.}
\item{plotfn}{If \code{plotfn = true}, then a plot of the noisy data and estimated signal are produced.}
\item{TI}{If TI = T, then the non-decimated transform is used. See the help pages for wd and wst for more on the non-decimated transform.}
\item{details}{If \code{details = FALSE} (the default), only the estimate of the underlying signal is returned. If \code{details = TRUE}, many other details are also returned.}
\item{policy}{Controls the type of thresholding done. Available policies are multiwavelet style (policy = "mws") and empirical Bayes (policy = "ebayes").}
\item{code}{Tells cthresh whether external C or NAG code is available to help with the calculations.}
\item{tol}{A tolerance parameter used in searching for prior parameters if the empirical Bayes policy is used.}
}
\details{
If a real-valued signal is decomposed using a complex-valued wavelet (like the Lina-Mayrand wavelets supplied by filter.select), then the wavelet coefficients are also complex-valued. Wavelet shrinkage can still be used to estimate the signal, by asking the question "which coefficients are small (and represent noise) and which are large (and represent signal)?" Two methods of determining which coefficients are small and which are large are proposed by Barber & Nason (2004). One is "multiwavelet style" thresholding (similar to that in Downie & Silverman (1998) where the coefficients are treated like the coefficients of a multiwavelet). Here, the "size" of the wavelet coefficient is determined as modulus of a standardised version of the coefficient. The standardisation is by the square root of the covariance matrix of the coefficient. A Bayesian method is to place a mixture prior on each coefficient. The prior has two components: a bivariate normal and a point mass at (0,0). The parameters are determined by an empirical Bayes argument and then the prior is updated by the data. 
}

\value{
Either a vector containing the estimated signal (if details = FALSE), or a list with the following components: 
\item{data}{The original data as supplied to cthresh.}
\item{data.wd}{The wavelet decomposition of the data.}
\item{thr.wd}{The thresholded version of data.wd.}
\item{estimate}{The estimate of the underlying signal.}
\item{Sigma}{The covariance matrices induced by the wavelet transform. See \code{make.dwwt} for more details.}
\item{sigsq}{The estimate of the variance of the noise which corrupted the data.}
\item{rule}{Which thresholding rule was used} 
\item{EBpars}{The empirical Bayes parameters found by the function find.parameters. Only present if the "ebayes" policy was used.}
\item{wavelet}{A list with components filter.number and family which, when supplied to \code{link{filter.select}}, determine the wavelet used to decompose the data.}} 

\note{
The estimates returned by cthresh have an imaginary component. In practice, this component is usually negligible. 
}
\section{RELEASE}{
Part of the CThresh addon to WaveThresh. Copyright Stuart Barber and Guy Nason 2004.}
\seealso{
\code{\link{filter.select}}, \code{\link{find.parameters}}, \code{\link{make.dwwt}}, \code{\link{test.dataCT}}, and the undocumented functions in CThresh. 

}
\examples{
#
# Make up some noisy data
#
y <- example.1()$y
ynoise <- y + rnorm(512, sd=0.1)
#
# Do complex-valued wavelet shrinkage with decimated wavelets
#
est1 <- cthresh(ynoise, TI=FALSE)
#
# Do complex-valued wavelet shrinkage with nondecimated wavelets
#
est2 <- cthresh(ynoise, TI=TRUE)
#
#
#
plot(1:512, y, lty=2, type="l")
lines(1:512, est1, col=2)
lines(1:512, est2, col=3)
}
\author{Stuart Barber}
\keyword{manip}
