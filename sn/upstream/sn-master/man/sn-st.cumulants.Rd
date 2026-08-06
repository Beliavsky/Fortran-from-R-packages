%  file sn/man/sn-st.cumulants.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013 Adelchi Azzalini
%---------------------
\name{sn-st.cumulants}
\alias{sn.cumulants}
\alias{st.cumulants}
\concept{cumulant}
\title{Cumulants of univariate skew-normal and skew-\eqn{t} distributions}

\description{Compute cumulants of univariate (extended) skew-normal and 
  skew-\eqn{t} distributions up to a given order.}

\usage{
  sn.cumulants(xi=0, omega=1, alpha=0, tau=0, dp=NULL, n=4)
  st.cumulants(xi=0, omega=1, alpha=0, nu=Inf, dp=NULL, n=4)
}
 
\arguments{
  \item{xi}{location parameters (numeric vector).}
  \item{omega}{scale parameters (numeric vector, positive).}
  \item{alpha}{slant parameters (numeric vector).}
  \item{tau}{hidden mean parameter (numeric scalar).}
  \item{nu}{degrees of freedom (numeric scalar, positive); the default value
    is \code{nu=Inf} which corresponds to the skew-normal distribution.}
  \item{dp}{a vector containing the appropriate set of parameters. 
    If \code{dp} is not \code{NULL}, the individual parameters must 
    not be supplied.}
  \item{n}{maximal order of the cumulants. For \code{st.cumulants} and
    for \code{sn.cumulants} with \code{tau!=0} (\acronym{ESN} distribution), 
    it cannot exceed 4.}
}

\section{Background}{
See Sections 2.1.4, 2.2.3 and 4.3.1 of the reference below}

\value{A vector of length \code{n} or a matrix with \code{n} columns, 
in case the input values are vectors.}

\references{
 Azzalini, A. with the collaboration of Capitanio, A. (2014). 
 \emph{The Skew-Normal and Related Families}. 
 Cambridge University Press, IMS Monographs series.
}

\author{Adelchi Azzalini}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{dsn}}, \code{\link{dsn}}}

\examples{
sn.cumulants(omega=2, alpha=c(0, 3, 5, 10), n=5)
sn.cumulants(dp=c(0, 3, -8), n=6)
st.cumulants(dp=c(0, 3, -8, 5), n=6) # only four of them are computed
st.cumulants(dp=c(0, 3, -8, 3))
}

\keyword{distribution}

