\name{find.parameters}
\alias{find.parameters}
\title{Find estimates of prior parameters}
\description{
Estimate the prior parameters for the complex empirical Bayes shrinkage procedure.
}
\usage{
find.parameters(data.wd, dwwt, j0, code, tol, Sigma)
}
\arguments{
        \item{data.wd}{Wavelet decomposition of the data being analysed.}
        \item{dwwt}{The diagonal elements of the matrix Wt(W). See \code{\link{make.dwwt}} for details.}
        \item{j0}{Primary resolution level, as discussed in the help for threshold.wd}
        \item{code}{Tells the function whether to use NAG code for the search (code="NAG"), R/S-plus for the search with C code to evaluate the likelihood (code="C"), or R/S-plus code for all calculations (code="R" or code="S"). Setting code="NAG" is strongly recommended.}
        \item{tol}{A tolerance parameter which bounds the mixing weight away from zero and one and the correlation between real and imaginary parts of the prior away from plus or minus one.}
        \item{Sigma}{The covariance matrix of the wavelet coefficients of white noise.}
}
\details{
The complex empirical Bayes (CEB) shrinkage procedure described by Barber & Nason (2004) places independent mixture priors on each complex-valued wavelet coefficient. This routine finds marginal maximum likelihood estimates of the prior parameters. If the NAG library is available, routine E04JYF is used otherwise the search is done using optimize (in R) or nlminb (in S-plus). In the latter case, the likelihood values should be computed externally using the C code supplied as part of the CThresh package - although a pure R / S-plus version is available, it is very slow. This function will not usually be called directly by the user, but is called from within cthresh.
}

\value{
A list with the following components:

        \item{pars}{Estimates of the prior parameters. Each row of this matrix contains the following parameter estimates for one level of the transform: mixing weight; variance of the real part of the wavelet coefficients; covariance between the real and imaginary parts; variance of the imaginary part of the wavelet coefficients. Note that for levels below the primary resolution, this search is not done and the matrix is full of zeros.}
        \item{Sigma}{The covariance matrix as supplied to the function.}
}

\note{
There may be warning messages from the NAG routine E04JYF. If the indicator variable IFAIL is equal to 5, 6, 7, or 8, then a solution has been found but there is doubt over the convergence. For IFAIL = 5, it is likely that the correct solution has been found, while IFAIL = 8 means that you should have little confidence in the parameter estimates. For more details, see the NAG software documentation available online at
\code{http://www.nag.co.uk/numeric/fl/manual19/pdf/E04/e04jyf_fl19.pdf}
}
\section{RELEASE}{
Part of the CThresh addon to WaveThresh. Copyright Stuart Barber and Guy Nason 2004.
}
\seealso{
\code{\link{cthresh}}
}
\author{Stuart Barber}
\keyword{manip}

