\name{sjpi}
\alias{sjpi}
\title{
Sheather-Jones Plug-in bandwidth criterion.
}
\usage{
sjpi(x, a)
}
\description{
  Given a dataset and set of pilot bandwidths, this function
  computes a bandwidth via the plug-in method, and the assumed
  `pilot' relationship of Sheather and Jones (1991).
  The S-J method chooses the bandwidth at which the two intersect.

  The purpose of this function is to demonstrate the sensitivity of
  plug-in methods to pilot bandwidths and assumptions. This function
  does not provide a reliable method of bandwidth selection.
}
\arguments{
  \item{x}{data vector}
  \item{a}{vector of pilot bandwidths}
}
\value{
  A matrix with four columns; the number of rows equals the length of \code{a}.
  The first column is the plug-in selected bandwidth. The second column
  is the pilot bandwidths \code{a}. The third column is the pilot bandwidth
  according to the assumed relationship of Sheather and Jones. The fourth
  column is an intermediate calculation.
}
\examples{
# Fig 10.2 (S-J parts) from Loader (1999).
data(geyser, package="locfit")
gf <- 2.5
a <- seq(0.05, 0.7, length=100)
z <- sjpi(geyser, a)

# the plug-in curve. Multiplying by gf=2.5 corresponds to Locfit's standard
# scaling for the Gaussian kernel.
plot(gf*z[, 2], gf*z[, 1], type = "l", xlab = "Pilot Bandwidth k", ylab
     = "Bandwidth h")

# Add the assumed curve.
lines(gf * z[, 3], gf * z[, 1], lty = 2)
legend(gf*0.05, gf*0.4, lty = 1:2, legend = c("Plug-in", "SJ assumed"))
}
\references{
  Sheather, S. J. and Jones, M. C. (1991). A reliable data-based bandwidth
  selection method for kernel density estimation. JRSS-B 53, 683-690.
}
\seealso{
  \code{\link{locfit}},
  \code{\link{locfit.raw}},
  \code{\link{lcvplot}}
}

\keyword{htest}
% Converted by Sd2Rd version 0.2-a5.
