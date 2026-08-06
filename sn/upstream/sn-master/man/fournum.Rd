%  file sn/man/fournum.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013-2020 Adelchi Azzalini
%---------------------
\name{fournum}
\encoding{UTF-8}
\alias{fournum}
\concept{quantile}
\concept{median}
\concept{variability}
\concept{asymmetry}
\concept{skewness}
\concept{kurtosis}
\title{Four-number summary of a numeric vector}
\description{
Returns a quantile-based four-number summary of the input data
}
\usage{
fournum(x, na.rm = TRUE, ...)
}
\arguments{
  \item{x}{a numeric vector,  maybe including \code{NA}s and \code{+/-Inf}'s.
     At least 8 not-\code{NA} values are required.
     It works with objects which can be coerced to vector.}
  \item{na.rm}{	 logical; 
     if \code{TRUE}, all \code{NA} and \code{NaN}s are dropped, 
     before the statistics are computed.}
 
  \item{\dots}{optional arguments passed to \code{\link[stats]{quantile}}}
}

\details{Function \code{quantile} is used to compute 7 octiles of \code{x}, 
  that is, quantiles of level \code{(1:7)/8}, denoted \code{oct[1:7]}, 
  and derive four summary quantities: 
  \enumerate{
  \item the median, which corresponds to \code{oct[4]}, 
  \item the \sQuote{(coefficient of) quartile deviation} 
    or semi-interquantile range: \code{(oct[6] - oct[2])/2};
  \item the Galton-Bowley measure of asymmetry, that is, skewness:
     \code{(oct[6] - 2 * oct[4] + oct[2])/(oct[6] - oct[2])};
  \item the Moors measure of kurtosis: 
       \code{(oct[7] - oct[5] + oct[3] - oct[1])/(oct[6] - oct[2])}
  }
  The term \sQuote{coefficient of quartile deviation} is adopted from the
  Encyclopedia of Statistical Sciences; see the references below.
  What is called Galton-Bowley measure here is often named 
  \sQuote{Bowley's measure}, but some sources attribute it to Francis Galton. 
  For the Moors measure, see Moors (1988).
}

\value{a vector of length four containing the median, the quartile deviation, 
  the Galton-Bowley measure and the Moors measure. However,
  if \code{x} does not contain at least 8 values (after removing \code{NA}s),
  \code{rep(NA,4)} is returned.}
  

\references{

\sQuote{Quartile deviation, coefficient of}, in:
  \emph{Encyclopedia of Statistical Sciences}, 2nd edition (2006). 
  Editors: Samuel Kotz (Editor-in-Chief), Campbell B. Read, N. Balakrishnan, 
  Brani Vidakovic. Volume 10, p.6743.
  
\sQuote{Skewness, Bowleys's measures of}, in:
  \emph{Encyclopedia of Statistical Sciences}, 2nd edition (2006). 
  Editors: Samuel Kotz (Editor-in-Chief), Campbell B. Read, N. Balakrishnan, 
  Brani Vidakovic.   Volume 12, p.7771-7773.


Moors, J.J.A. (1988). 
  A quantile alternative for kurtosis.
  \emph{Source: Journal of the Royal Statistical Society. 
  Series D (The Statistician)}, Vol. 37, pp. 25-32

}
\author{Adelchi Azzalini}

\note{
 Computation of octiles makes real sense only if \code{length(x)} is
 substantially larger than 8. }
  
\seealso{\code{\link[stats]{quantile}}, \code{\link[stats]{fivenum}},
   \code{\link[stats]{IQR}}  }

\examples{
fournum(datasets::rivers)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{univar}
\keyword{nonparametric} 
\keyword{robust} 
