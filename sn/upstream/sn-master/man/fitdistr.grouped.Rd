%  file sn/man/fitdistr.grouped.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2022 Adelchi Azzalini
%
\name{fitdistr.grouped}
\alias{fitdistr.grouped}
\title{
Maximum-likelihood fitting of a univariate distribution from grouped data
}
\description{
Maximum-likelihood fitting of a univariate distribution when the data are
represented by a set of frequencies pertaining to given set of contiguous 
intervals.
}
\usage{
fitdistr.grouped(breaks, counts, family, weights, trace = FALSE, wpar = NULL)
}
\arguments{
  \item{breaks}{A numeric vector of strictly increasing values which identify 
     a set of contiguous intervals on the real line.
     See \sQuote{Details} for additional information.}
     
  \item{counts}{A vector of non-negative integers representing the number 
     of observations falling in the intervals specified by \code{breaks};
     it is then required that \code{length(counts)+1=length(breaks)}. }
     
  \item{family}{A character string specifying the parametric family of 
     distributions to be used for fitted. 
     Admissible names are: \code{"normal"}, \code{"logistic"}, \code{"t"}, 
     \code{"Cauchy"}, \code{"SN"}, \code{"ST"}, \code{"SC"}, \code{"gamma"}, 
     \code{"Weibull"};
     the names \code{"gaussian"} and \code{"Gaussian"} are also allowed, 
     and are converted to \code{"normal"}.}
     
  \item{weights}{An alias for \code{counts}, allowed for analogy with the 
    \code{selm} function.}

  \item{trace}{A logical value which indicates whether intermediate evaluations 
    of the optimization process are printed (default: \code{FALSE}).}
    
  \item{wpar}{
    An optional vector with initial values of the \sQuote{working parameters}
    for starting the maximization of the log-likelihood function; 
    see \sQuote{Details} for their description. }
  }

\details{
The original motivation of this function was fitting a univariate \acronym{SN}, 
\acronym{ST} or \acronym{SC} distribution from grouped data;
its scope was later extended to include some other continuous distributions. 
The adopted name of the function reflects the broad similarity of its purpose 
with the one of \code{\link[MASS]{fitdistr}}, but there are substantial 
differences in the actual working of the two functions.

The parameter set of a given \code{family} is the same as appearing in the
corresponding \code{d<basename>} function, with the exception of the \code{"t"} 
distribution, for which a location and a scale parameter are included, 
besides \code{df}.
   
The range of \code{breaks} does not need to span the whole support of the
chosen \code{family} of distributions, that is, 
\code{(0, Inf)} for \code{"Weibull"} and \code{"gamma"} families,
\code{(-Inf, Inf)} for the other families.
In fact, for the purpose of post-fitting plotting, an infinite 
\code{range(breaks)} represents a complication, requiring an \emph{ad hoc} 
handling;  so it is sensible to avoid it.
However, at the maximum-likelihood fitting stage, the full support of
the  probability distribution is considered, with the following implications. 
If \code{max(breaks)=xR}, say, and \code{xR<Inf},  then an additional 
interval \code{(xR, Inf)} is introduced, with value \code{counts=0} assigned.
A similar action is taken at the lower end: if \code{min(breaks)=xL} is
larger than the infimum of the support of the distribution,
an extra  0-\code{counts} interval is introduced as \code{(0, xL)} 
or \code{(-Inf, xL)}, depending on the support of the \code{family}.

Maximum likelihood fitting is obtained by maximizing the pertaining 
multinomial log-likelihood using the \code{\link[stats]{optim}} function 
with  \code{method="Nelder-Mead"}. For numerical convenience, the numerical 
search is performed using \sQuote{working parameters} in place of the original
ones, with reverse conversion at the end. The working parameters coincide
with the original distribution parameters when they have unbounded range,
while they are log-transformed in case of intrinsically positive parameters.
This transformation applies to the parameters of the positive-valued
distributions ("gamma" and "Weibull"), all scale parameters and \code{df}
of the \code{"t"} distribution.
}

\value{
An object of class \code{fitdistr.grouped}, whose components are described in 
\code{\link{fitdistr.grouped-class}}.
}

%%  \references{multinomial likelihood?}

\author{Adelchi Azzalini }
%%\note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{For methods pertaining to this class of objects, see
\code{\link{fitdistr.grouped-class}} and
\code{\link{plot.fitdistr.grouped}}; see also 
\code{\link{dsn}}, \code{\link{dst}}, \code{\link{dsc}},
\code{\link[stats]{Distributions}}, \code{\link[stats]{dmultinom}};
see also \code{\link{selm}} for ungrouped data fitting and an example
elaborated on below.
}

\examples{
data(barolo)
attach(barolo)
A75 <- (reseller=="A" & volume==75)
logPrice <- log(price[A75], 10) # original data fitting, as in documentation of 'selm' 
detach(barolo)
#  now construct grouped data  
breaks<- seq(1, 3, by=0.25)
f <- cut(logPrice, breaks = breaks)
counts <- tabulate(f, length(levels(f))) 
# grouped data fitting
logPrice.grouped <- fitdistr.grouped(breaks, counts, family='ST') 
summary(logPrice.grouped) # compare this fit with the ungrouped data fitting 
}

\keyword{distribution}
\keyword{univar}
\concept{grouped data}
\concept{multinomial distribution} 