%  file sn/man/makeSUNdistr.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013-2021 Adelchi Azzalini
%---------------------
\name{makeSUNdistr}
\docType{package}
\encoding{UTF-8}

\alias{makeSUNdistr}

\title{Build an object representing a SUN distribution}

\description{
Build an object which identifies a Unified Skew-Normal distribution 
(\acronym{SUN}) within this parametric family.
The \acronym{SUN} family is essentially equivalent to some other 
parametric families examined in the literature, notably the Closed Skew-Normal.
}

\usage{makeSUNdistr(dp, name, compNames, HcompNames, drop = TRUE)}

\arguments{
   \item{dp}{a list of parameters as described at \code{\link{SUNdistr-base}}.}
  \item{name}{an optional character string with the name of the distribution. 
      If missing, one is created.}
  \item{compNames}{an optional vector of character strings with the names of 
      the component variables; its length must be equal to the dimensionality 
      \code{d} of the distribution being generated.
      If missing, the components are named  \code{"V1"}, \code{"V2"}, \dots}
  \item{HcompNames}{an optional vector of character strings with the names of 
      the hidden component variables; its length must be equal to the 
      dimensionality component \code{m} described in the \sQuote{Details}.
      If missing, the components are named \code{"H1"}, \code{"H2"}, \dots}
  \item{drop}{a logical value (default: \code{TRUE}) relevant only in the
      case \code{m=1}. When both \code{m=1} and \code{drop=TRUE}, the 
      returned object is of class either \code{SECdistrUv} or \code{SECdistrMv},
      depending on the value of \code{d}, and family \code{"SN"} or 
      \code{"ESN"}, depending on the \code{dp} ingredients.}
}

\details{
The argument  \code{dp} is a list, whose components are described at 
\code{\link{SUNdistr-base}}; see especially the \sQuote{Details} there.
In this respect, there is no difference between the univariate and the
univariate case, differently from the similar command 
\code{\link{makeSECdistr}}.

If the arguments \code{name}, \code{compNames} and \code{HcompNames}
are missing, they are composed from the supplied arguments.

A \code{SUNdistr-class} object operates according to the S4 protocol.
}

\value{An object of \code{\link{SUNdistr-class}}}

\author{Adelchi Azzalini}

\note{
The present structure and user interface of this function, and of other ones 
related to the \acronym{SUN} distribution, must be considered experimental, 
and they might possibly change in the future.} 

\seealso{
 Basic information on the SUN distribution \code{\link{SUNdistr-base}},
 the description of the class \code{\link{SUNdistr-class}}, 

Related methods:
\code{\link{show.SUNdistr}} for displaying the object constituents,
\code{\link{plot.SUNdistr}} for plotting,
\code{\link{mean.SUNdistr}} for the mean value,
\code{\link{vcov.SUNdistr}} for the variance matrix,
\code{\link{summary.SUNdistr}} for various summary quantities 

Functions \code{\link{SUNdistr-op}}  manipulate objects created by this 
function, producing new \code{SUNdistr-class} objects  
}

\examples{
xi <- c(1, 0, -1)
Omega <- matrix(c(2,1,1, 1,3,1, 1,1,4), 3, 3)
Delta <- matrix(c(0.72,0.20, 0.51,0.42, 0.88, 0.94), 3, 2, byrow=TRUE)
Gamma <- matrix(c(1, 0.8, 0.8, 1), 2, 2)
dp3 <- list(xi=xi, Omega=Omega, Delta=Delta, tau=c(-0.5, 0), Gamma=Gamma)
sun3 <- makeSUNdistr(dp=dp3, name="SUN3", compNames=c("x", "w", "z"))
show(sun3)
}

\keyword{distribution}
\keyword{multivariate} 
\concept{SUN distribution}
\concept{Unified Skew-Normal distribution}
\concept{CSN distribution}
\concept{Closed Skew-Normal distribution}
