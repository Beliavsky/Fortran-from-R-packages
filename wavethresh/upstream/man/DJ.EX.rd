\name{DJ.EX}
\alias{DJ.EX}
\title{Produce Donoho and Johnstone test functions}
\usage{
DJ.EX(n=1024, signal=7, rsnr=7, noisy=FALSE, plotfn=FALSE)
}
\arguments{
\item{n}{Number of samples of the function required.}
\item{signal}{A factor that multiples the function values.}
\item{rsnr}{If Gaussian noise is to be added to the functions
	then this argument specifies the root signal to noise ratio.}
\item{noisy}{If TRUE then Gaussian noise is added to the signal so that the root signal to noise ratio is \code{rsnr}.
If FALSE then just the signals are returned.}
\item{plotfn}{If TRUE then a plot is produced. If FALSE no plot is produced.}
}
\description{
Function to produce the blocks, bumps, Doppler and heavisine functions described by Donoho and Johnstone (1994).

}
\details{
The Donoho and Johnstone test functions were designed to reproduce various features to be found in real world signals such as jump discontinuities (blocks), spikes (the NMR-like bumps), varying frequency behaviour (Doppler) and jumps/spikes in smooth signals (heavisine).
These functions are most often used for testing wavelet shrinkage methods and comparing them to other nonparametric regression techniques.
(Donoho, D.L. and Johnstone, I.M. (1994), Ideal spatial adaptation by wavelet
shrinkage. \emph{Biometrika}, \bold{81}, 425--455).

Another version of the Doppler function can be found in the standalone
\code{\link{doppler}} function.

Another function for this purpose is the Piecewise Polynomial created in Nason and Silverman (1994) an encapsulated in WaveThresh by \code{\link{example.1}}
(Nason, G.P. and Silverman, B.W. (1994) The discrete wavelet transform in
S, \emph{J. Comput. Graph. Statist.}, \bold{3}, 163--191.

\emph{NOTE: This function might not give exactly the same function values as the equivalent function in WaveLab}
}
\value{
A list with four components: blocks, bumps, heavi and doppler containing the sampled signal values for the four types of Donoho and Johnstone test functions. Each of these are deemed to be sampled on an equally spaced grid from 0 to 1.
}
\seealso{\code{\link{doppler}},\code{\link{example.1}}, \code{\link{threshold}}, \code{\link{wd}}}
\examples{
#
# Show a picture of the four test functions with the default args
#
\dontrun{DJ.EX(plotfn=TRUE)}
}
\author{Theofanis Sapatinas}
\keyword{nonparametric}
