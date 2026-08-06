%  file sn/man/summary.SECdistr.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013-2015 Adelchi Azzalini
%---------------------
\name{summary.SECdistr}
\alias{summary.SECdistr}
\alias{summary.SECdistrUv}
\alias{summary.SECdistrMv}
\alias{summary,SECdistrUv-method}
\alias{summary,SECdistrMv-method}

\title{Summary of a \acronym{SEC} distribution object}

\description{Produce a summary of an object of class either
\code{"SECdistrUv"} or \code{"SECdistrMv"}, which refer to a univariate or a
multivariate \acronym{SEC} distribution, respectively. Both types of
objects can be produced by \code{makeSECditr}. }

\usage{
\S4method{summary}{SECdistrUv}(object, cp.type = "auto", probs)

\S4method{summary}{SECdistrMv}(object, cp.type = "auto")
}

\arguments{
  \item{object}{an object of class \code{"SECdistrUv"} or \code{"SECdistrMv"}.}

  \item{cp.type}{a character string to select the required variance of
    \acronym{CP} parameterization; possible values are \code{"proper"},
    \code{"pseudo"}, \code{"auto"} (default). For a description of these
    codes, see \code{\link{dp2cp}}.}
  
 \item{probs}{in the univariate case, a vector of probabilities for which 
    the corresponding quantiles are required. If missing, the vector
     \code{c(0.05, 0.25, 0.50, 0.75, 0.95)} is used.}
}


\details{For a description of  the \acronym{DP},  \acronym{CP} and    
pseudo-\acronym{CP}  parameter sets included in the returned object, 
see \code{\link{dp2cp}}.

The \code{aux} slot of the returned object includes other summary quantities,
as described next. 
In the univariate case, the reported quantile-based measures of skewness and
kurtosis refer to the  Bowley and Moors measures, respectively;
see Groeneveld (2006) and Moors (1988) for their specifications.
In the multivariate case, the Mardia's measures of skewness and kurtosis
are computed from the expressions given on p.153 and p.178 of
Azzalini and Capitanio (2014). 
Note that the Mardia measures are computed in their \eqn{\gamma}-form, not the 
\eqn{\beta}-form; in practice, these only differ by a \eqn{d(d+2)} shift 
in the measure of kurtosis, if \eqn{d} denotes the dimension.

In the univariate case, \code{delta} is a simple transformation of the 
slant parameter \code{alpha}; it takes values in \eqn{(-1, 1)}. 
In the multivariate case, \code{delta} is a vector with components
of similar type; they correspond to the matching terms of the univariate 
components. The \code{alpha*} and \code{delta*} coefficients are univariate
comprehensive summary quantities of slant; see pp.132-3 of 
Azzalini and Capitanio (2014) for their expressions. These quantities
play an important role in \acronym{SEC} distributions; for instance,
the Mardia's measures of multivariare skewness and kurtosis depend
on the vector of slant parameters only via \code{delta*} or, equivalently,
via \code{alpha*}.

The mode, which is unique for all these distributions, is computed by a 
numerical line search between the \acronym{DP} location and the \acronym{CP}
location (or the pseudo-\acronym{DP} location, when the latter does
exists). This line search is univariate also in the multivariate case, 
using Propositions 5.14 and 6.2 of Azzalini and Capitanio (2014);
see also Problem 5.14.

The examples below illustrate how extract various components from \code{aux} 
and other slots of the returned object. 
}

 

\value{A list with the following components:
  \item{family}{name of the family within the \acronym{SEC} class, character}
  \item{dp}{\acronym{DP} parameters, a list or a vector}
  \item{name}{the name of the distribution, character string}
  \item{compNames}{in the multivariate case the names of the components,
     a character vector}
  \item{cp}{\acronym{CP} parameters, a list or a vector}
  \item{cp.type}{the name of the selected variant of the \acronym{CP} set}
  \item{aux}{a list with auxiliary ingredients (mode, coefficients of 
       skewness and kurtosis, in the parametric and non-parametric variants, 
       and more); see Section \sQuote{Details} for more information.}
       
The list items \code{dp} and \code{cp} are vectors if \code{class(object)} is  
\code{SECdistrUv} (univariate distribution); they are lists if 
\code{class(object)} is \code{SECdistrMv} (multivariate distribution). 

}
 
\author{Adelchi Azzalini}

%\note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{makeSECdistr}} for building a \acronym{SEC} distribution

 \code{\link{extractSECdistr}} for extracting a \acronym{SEC}
 distribution from a \code{\link{selm}} fit 
 
 methods \code{\link[base]{mean}} and \code{\link[sn]{sd}} 
 for computing the mean  and the standard deviation of 
 \code{\link{SECdistrUv-class}} objects,
 methods \code{\link[base]{mean}} and  \code{\link[stats]{vcov}} 
 for computing the mean vector and the variance matrix of 
 \code{\link{SECdistrMv-class}} objects
 
 \code{\link[sn]{modeSECdistr}} for computing the mode directly 
 }
 
 
\references{

Azzalini, A. with the collaboration of Capitanio, A. (2014). 
 \emph{The Skew-Normal and Related Families}. 
 Cambridge University Press, IMS Monographs series.

Moors, I. J. A. (1988). A quantile alternative for kurtosis. 
\emph{The Statistician} \bold{37}, 25-32.

Groeneveld, R. A. (2006).  Skewness, Bowley's measures of. 
In volume \bold{12}, 7771-3, of \emph{Encyclopedia of Statistical Sciences}, 
2nd edition, edited by Kotz et al. Wiley, New York.
}


\examples{
f3 <- makeSECdistr(dp=c(3,2,5), family="SC")
summary(f3)
s <- summary(f3, probs=(1:9)/10)
print(slotNames(s)) 
print(names(slot(s,"aux")))   # the components of the 'aux' slot
slot(s, "aux")$mode           # the same of modeSECdistr(object=f3)
slot(s, "aux")$q.measures     # quantile-based measures of skewness and kurtosis
#
dp3 <- list(xi=1:3, Omega=toeplitz(1/(1:3)), alpha=c(-3, 8, 5), nu=6)
st3 <- makeSECdistr(dp=dp3, family="ST", name="ST3", compNames=c("U", "V", "W"))
s <- summary(st3)
dp <- slot(s, "dp")           # the same of slot(st3, "dp")
slot(s, "cp")$var.cov         # the same of vcov(st3)
slot(s, "aux")$delta.star     # comprehensive coefficient of shape
slot(s, "aux")$mardia         # Mardia's measures of skewness and kurtosis
#
dp2 <- list(xi=rep(0,2), Omega=matrix(c(2,2,2,4),2,2), alpha=c(3,-5), tau=-1)
esn2 <- makeSECdistr(dp=dp2, family="ESN", name="ESN-2d")
summary(esn2)
}
\keyword{multivariate}
\keyword{distribution} 
