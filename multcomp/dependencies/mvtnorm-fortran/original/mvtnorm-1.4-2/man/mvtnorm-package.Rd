\name{mvtnorm-package}
\alias{mvtnorm-package}
\alias{mvtnorm}
\docType{package}
\title{
\packageTitle{mvtnorm}
}
\description{
\packageDescription{mvtnorm}
}
\details{
Package \pkg{mvtnorm} provides functionality for dealing with multivariate
normal and t-distributions. The package interfaces \proglang{FORTRAN} and
\code{C} code for evaluating multivariate normal probabilities written by
Alan Genz and Tetsuhisa Miwa. Functions \code{\link{pmvnorm}},
\code{\link{pmvt}}, \code{\link{qmvnorm}}, and \code{\link{qmvt}} return
normal and t probabilities or corresponding quantiles computed by these
original implementations. Users interested in the computation of such
probabilities or quantiles, for example for multiple testing purposes,
should use this functionality.

When the multivariate normal log-likelihood function, defined by the
log-probability in the discrete or interval-censored case or by the
log-density for exact real observations, or a mix thereof, shall be
computed, functions \code{\link{lpmvnorm}}, \code{\link{ldmvnorm}}, and
\code{\link{ldpmvnorm}} are better suited.  They rely on an independent
implementation of Genz' algorithm (for log-probabilities), can be customised
(different quasi-Monte Carlo schemes), and are a bit faster.  Most
importantly, the corresponding score functions are available through
functions \code{\link{slpmvnorm}}, \code{\link{sldmvnorm}}, or
\code{\link{sldpmvnorm}}, which help to speed-up parameter estimation
considerably. Users interested in this functionality should
consult the \code{lmvnorm_src} package vignette.

}
\seealso{\code{vignette("lmvnorm_src", package = "mvtnorm")}}
\keyword{ package }