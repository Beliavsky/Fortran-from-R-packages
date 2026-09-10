\name{doppler}
\alias{doppler}
\title{Evaluate the Donoho and Johnstone Doppler signal.
}
\usage{
doppler(t)
}
\arguments{
\item{t}{The domain of the Doppler function (where you wish to evaluate this
	Doppler function}
}
\description{
This function evaluates and returns the Doppler signal from Donoho and Johnstone, (1994).
}
\details{
This function evaluates and returns the Doppler signal from Donoho and Johnstone, (1994).
(Donoho, D.L. and Johnstone, I.M. (1994), Ideal spatial adaptation by wavelet
shrinkage. \emph{Biometrika}, \bold{81}, 425--455).

Another version of this function can be found in \code{\link{DJ.EX}}.
}
\value{
A vector of the same length as the input vector containing the Doppler signal
at \code{t}
}
\seealso{\code{\link{DJ.EX}}}
\examples{
#
# Evalute the Doppler signal at 100 arbitrarily spaced points.
#
tt <- sort(runif(100))
dopp <- doppler(tt)
\dontrun{plot(tt, dopp, type="l")}
}

\author{Arne Kovac}
\keyword{nonparametric}
