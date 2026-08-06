%  file sn/man/convertCSN2SUNpar.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013-2021 Adelchi Azzalini
%---------------------
\name{convertCSN2SUNpar}
\encoding{UTF-8}

\alias{convertCSN2SUNpar}

\title{Conversion of \acronym{CSN} parameters to \acronym{SUN} parameters} 

\description{The parameter set of a Closed Skew-Normal (\acronym{CSN})
distribution is  converted into the parameter set of the equivalent 
Unified Skew-Normal (\acronym{SUN}) distribution.}

\usage{convertCSN2SUNpar(mu, Sigma, D, nu, Delta)}

\arguments{
  \item{mu}{a numeric vector of length \code{p}, say.}
  \item{Sigma}{a positive definite variance matrix of size \code{c(p,p)}.}
  \item{D}{an arbitrary numeric matrix of size say \code{c(q, p)}, say.}
  \item{nu}{a numeric vector of length \code{q}.}
  \item{Delta}{a positive definite variance matrix of size \code{c(q,q)}.}
}

\details{The arguments of the function match the parameters 
\eqn{(\mu, \Sigma, D, \nu, \Delta)} of the \acronym{CSN} distribution
presented by González-Farías \emph{et alii} (2004a, 2004b). 
These parameters are converted into those of the equivalent 
\acronym{SUN} distribution, which is unique. The converse operation, that is,
mapping parameters from the \acronym{SUN} to the \acronym{CSN} family, 
is not handled here. Its solution  would be non-unique, 
because the \acronym{CSN} family is over-parameterized.

Note that, having retained the exact notation of the above-quoted papers,
there is a \code{Delta} argument which must not be confused with one of the
arguments for the \acronym{SUN} distribution  in \code{\link{SUNdistr-base}}. 
The coincidence of these names is entirely accidental.

The \acronym{CSN} parameters must only satisfy the requirements that 
\eqn{\Sigma} and \eqn{\Delta} are symmetric positive definite matrices. 
Since these conditions are somewhat simpler to check than those for the 
\acronym{SUN} parameters, as indicated in \code{\link{SUNdistr-base}}, 
this function may provide a simple option for the specification of a
\acronym{CSN/SUN} distribution.

The parameter list \code{dp} produced by this function can be used as an
input for the functions in \code{\link{SUNdistr-base}} or for
\code{\link{makeSUNdistr}}.
}

\value{a list representing the \code{dp} parameter set of the
       corresponding \acronym{SUN} distribution}

\references{
González-Farías, G., Domínguez-Molina, J. A., & Gupta, A. K. (2004a). 
Additive properties of skew normal random vectors.
\emph{J. Statist. Plann. Inference} \bold{126}, 521-534.

% González-Farías, G., Domínguez-Molina, J. A., & Gupta, A. K. (2004b). 
González-Farías, G., Domínguez-Molina, J. A., & Gupta, A. K. (2004b). 
The closed skew-normal distribution. 
In M. G. Genton (Ed.), \emph{Skew-elliptical Distributions and Their 
Applications: a Journey Beyond Normality}, Chapter 2, (pp. 25–42). 
Chapman & Hall/\acronym{CRC}. 
}
\author{Adelchi Azzalini}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{SUNdistr-base}}, \code{\link{makeSUNdistr}} }

\examples{
p <- 3
q <- 2
mu <- 1:p
Sigma <- toeplitz(1/(1:p))
D <- matrix(sqrt(1:(p*q)), q, p)
nu <- 1/(1:q)
Delta <- diag(q) + outer(rep(1,q), rep(1,q))
dp <- convertCSN2SUNpar(mu, Sigma, D, nu, Delta)
}

\keyword{distribution}
\keyword{multivariate} 
\concept{SUN distribution}
\concept{Unified Skew-Normal distribution}
\concept{CSN distribution}
\concept{Closed Skew-Normal distribution}
 
 
