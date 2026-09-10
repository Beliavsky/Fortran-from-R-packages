\name{threshold.mwd}
\alias{threshold.mwd}
\title{Use threshold on an mwd object. }
\description{
Applies hard or soft thresholding to multiple wavelet decomposition object mwd.object.  
}
\usage{
\method{threshold}{mwd}(mwd, levels = 3:(nlevelsWT(mwd) - 1), type = "hard",
    policy = "universal", boundary = FALSE, verbose = FALSE,
    return.threshold = FALSE, threshold = 0, covtol = 1e-09, 
    robust = TRUE, return.chisq = FALSE,
    bivariate = TRUE, \dots)
}
\arguments{
\item{mwd}{The multiple wavelet decomposition object that you wish to threshold.}
\item{levels}{a vector of integers which determines which scale levels are thresholded in the decomposition. Each integer in the vector must refer to a valid level in the \code{\link{mwd}} object supplied. This is usually any integer from 0 to \code{\link{nlevelsWT}}(wd)-1 inclusive. Only the levels in this vector contribute to the computation of the threshold and its application. }
\item{type}{determines the type of thresholding this can be "\code{hard}" or "\code{soft}".}
\item{policy}{selects the technique by which the threshold value is selected. Each policy corresponds to a method in the literature. At present the different policies are "\code{universal}", "\code{manual}", "\code{single}". The policies are described in detail below. }
\item{boundary}{If this argument is \code{TRUE} then the boundary bookeeping values are included for thresholding, otherwise they are not. }
\item{verbose}{if \code{TRUE} then the function prints out informative messages as it progresses. }
\item{return.threshold}{If this option is \code{TRUE} then the actual \emph{value} of the threshold is returned. If this option is FALSE then a thresholded version of the input is returned.}
\item{threshold}{This argument conveys the user supplied threshold. If the \code{policy="manual"} then \code{value} is the actual threshold value. Any other \code{policy} means that the \code{threshold} value is ignored.}
\item{covtol}{The tolerance for what constitutes a singular variance matrix. If smallest eigenvalue of the estimated variance matrix is less than \code{covtol} then it is assumed to be singular and no thresholding is done at that level. Note: do not confuse \code{covtol} with \code{cvtol} an argument in \code{\link{threshold.wd}}.} 
\item{robust}{If TRUE the variance matrix at each level is estimated using a robust method (mad) otherwise it is estimated using var().}
\item{return.chisq}{If TRUE the vector of values to be thresholded is returned. These values are a quadratic form of each coefficient vector, and under normal assumptions the noise component will have a chi-squared distribution (see Downie and Silverman 1996). }
\item{bivariate}{this line is in construction}
\item{\dots}{any other arguments}
}
\details{
Thresholding modifies the coefficients within a \code{\link{mwd.object}}. The modification can be performed either with a "hard" or "soft" thresholding selected by the type argument. 

Unless policy="single", the following method is applied. The columns of
\code{mwd$D} are taken as coefficient vectors \eqn{D_{j,k}}.
From these
\eqn{\chi^2_{j,k}=D_{j,k} \cdot V_j^{-1}}.
\eqn{D_{j,k}} is computed, where
\eqn{V_j^{-1}} is the inverse of the estimated variance of the coefficient vectors in that level (j).
\eqn{\chi^2_{j,k}} is a positive scalar which is to be thresholded in a
similar manner to univariate hard or soft thresholding.
To obtain the new values of \eqn{D_{j,k}} shrink the vector by the
same proportion as was the corresponding \eqn{\chi^2_{j,k}} term. 
i
}
\value{
An object of class \code{\link{mwd}}. This object contains the thresholded wavelet coefficients. Note that if the \code{return.threshold} option is set to TRUE then the threshold values will be returned, or if \code{return.chisq} the vector of values to be thresholded will be returned, rather than the thresholded object.} 

\note{

POLICIES 
\describe{
\item{single}{If \code{policy="single"} then univariate thresholding is applied to each element of D as in (Strela et al 1999).} 
\item{universal}{The \code{universal} threshold is computed using 2log(n) (See Downie & Silverman 1996) where n is the number of coefficient vectors to be thresholded.}
\item{manual}{The "\code{manual}" policy is simple. You supply a \code{threshold} value to the threshold argument and hard or soft thresholding is performed using that value}
}
}
\section{RELEASE}{Version 3.9.6 (Although Copyright Tim Downie 1995-6). }
\seealso{
\code{\link{accessC.mwd}}, \code{\link{accessD.mwd}}, \code{\link{draw.mwd}}, \code{\link{mfirst.last}}, \code{\link{mfilter.select}}, \code{\link{mwd}}, \code{\link{mwd.object}}, \code{\link{mwr}}, \code{\link{plot.mwd}}, \code{\link{print.mwd}}, \code{\link{putC.mwd}}, \code{\link{putD.mwd}}, \code{\link{summary.mwd}}, \code{\link{wd}}, \code{\link{wr.mwd}}. 
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
##
# Plot it
#
\dontrun{ts.plot(ynoise)}
#
# Now take the discrete multiple wavelet transform
# N.b. I have no idea if the default wavelets here are appropriate for
# this particular examples. 
#
ynmwd <- mwd(ynoise)
\dontrun{plot(ynwd)}
# [1] 2.020681 2.020681 2.020681 2.020681 2.020681 2.020681 2.020681
#
# Now do thresholding. We'll use the default arguments.
#
ynmwdT <- threshold(ynmwd)
#
# And let's plot it
#
\dontrun{plot(ynmwdT)}
#
# Let us now see what the actual estimate looks like
#
ymwr <- wr(ynmwdT)
#
# Here's the estimate... 
#
\dontrun{ts.plot(ywr1)}
}
\keyword{smooth}
\keyword{nonlinear}
\author{Tim Downie}
