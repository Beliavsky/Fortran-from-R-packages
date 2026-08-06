%  file sn/man/dp2cp.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013 Adelchi Azzalini
%---------------------
\name{dp2cp}
\alias{dp2cp}
\alias{cp2dp}
\alias{dp2op}
\alias{op2dp}
\title{Conversion between parametrizations of a skew-elliptical distribution}

\description{
  Convert direct parameters (\acronym{DP}) to centred parameters 
 (\acronym{CP}) of a skew-elliptical distribution and \emph{vice versa}.}
  
\usage{
dp2cp(dp, family, object = NULL, cp.type = "proper", upto = NULL) 
cp2dp(cp, family)
dp2op(dp, family)
op2dp(op, family)
}
 
\arguments{

\item{dp}{a vector (in the univariate case) or a list (in the multivariate
  case) as described in \code{\link{makeSECdistr}}; see \sQuote{Background 
  and Details} for an extented form of usage.}

\item{cp}{a vector or a list, in agreement with \code{dp} as for type and
  dimension.}

\item{op}{a vector or a list, in agreement with \code{dp} as for type and
  dimension.}
 
\item{family}{a characther string with the family acronym, 
   as described in \code{\link{makeSECdistr}}, except that family 
  \code{"ESN"} is not implemented.}

\item{object}{optionally, an S4 object of class \code{SECdistrUv} or
  \code{SECdistrMv}, as produced by \code{\link{makeSECdistr}} 
  (default value: \code{NULL}). 
  If this argument is not \code{NULL}, then \code{family} and \code{dp} 
  must not be set.}

\item{cp.type}{character string, which has effect only if \code{family="ST"} 
   or \code{"SC"}, otherwise a warning message is generated.  Possible values 
   are \kbd{"proper", "pseudo", "auto"}, which correspond to the \acronym{CP} 
   parameter set, their `pseudo-\acronym{CP}' version and an automatic
   selection based on \code{nu>4}, where \code{nu} represents the degrees of
   freedom of the \acronym{ST} distribution.}

\item{upto}{numeric value (in \code{1:length(dp)}, default=\code{NULL}) to 
  select how many \acronym{CP} components are computed. 
  Default value \code{upto=NULL} is equivalent to \code{length(dp)}.}
 } 

\value{For \code{dp2cp}, a matching vector (in the univariate case) or a list 
   (in the multivariate case) of \code{cp} parameters. 
   For \code{cp2dp} and \code{op2dp}, a similar object of \code{dp} parameters,
   provided the set of input parameters is in the admissible region.
   For \code{dp2op}, a similar set of \code{op} parameters.}
		
\section{Background}{For a description of the \acronym{DP}
parameters, see Section \sQuote{Details} of \code{\link{makeSECdistr}}.  The
\acronym{CP} form of parameterization is cumulant-based. For a univariate
distribution, the \acronym{CP} components are the mean value (first cumulant), 
the standard deviation (square root of the 2nd cumulant), the coefficient of 
skewness (3rd standardized cumulant) and,  for the \acronym{ST}, 
the coefficient of excess kurtosis (4th standardized cumulant). 
For a multivariate distribution, there exists an extension based on the 
same logic; its components represent the
vector mean value, the variance matrix, the vector of marginal coefficients of
skewness and, only for the \acronym{ST}, the Mardia's coefficient of excess
kurtosis. The pseudo-\acronym{CP} variant provides an `approximate form' of
\acronym{CP} when not all required cumulants exist; however, this parameter set
is not uniquely invertible to \acronym{DP}. The names of pseudo-\acronym{CP}
components printed in summary output are composed by adding a \code{~} 
after the usual component name; for example, the first one is denoted
\code{mean~}.

Additional information is provided by Azzalini and Capitanio (2014).
Specifically, their Section 3.1.4 presents \acronym{CP} in the univariate 
\acronym{SN} case, Section 4.3.4 \acronym{CP} for the \acronym{ST} case and 
the `pseudo-\acronym{CP}' version. Section 5.2.3 presents the multivariate
extension for the \acronym{SN} distribution, Section 6.2.5 for the 
multivariate \acronym{ST} case. 
For a more detailed discussion, see Arellano-Valle & Azzalini (2013).

In the univariate case, the \acronym{OP} parameterization is very similar to 
\acronym{DP}, from which  it differs only for the components which regulate 
dispersion (or scatter) and slant.  
The real of relevance the \acronym{OP} parameterization appears in the 
multivariate case, where  the components of the slant parameter can be 
interpreted component-wise and remain unaffected if marginalization 
with respect to some other components  is performed. 
In the multivariate \acronym{SN} case, the components of \acronym{OP}, denoted
\eqn{\xi, \Psi, \lambda}, are associated to the expression of the density
function (5.30) of Azzalini & Capitanio (2014); see the discussion on p.131 
for more information. In the univariate case, the slant components of 
\acronym{DP} and \acronym{OP} coincide, that is, \eqn{\alpha=\lambda}.
The acronym \acronym{OP} stands for `original parameterization' since this is, 
up to a very minor difference, the parameterization adopted by Azzalini & 
Dalla Valle (1996).
}

\section{Details}{
While any choice of the components of \acronym{DP} or \acronym{OP} is
admissible, this is not true for \acronym{CP}. An implication is that a
call to \code{cp2dp} may fail with an error message \code{"non-admissible CP"}
for certain input values. The most extreme case is represented by the
\acronym{SC} family, for which  \acronym{CP} never exists; hence it makes
to sense to call \code{cp2dp} with \code{family="SC"}.

It is possible to call the functions with \code{dp} or \code{cp} having more
components than those expected for a given family as described above and in
\code{\link{makeSECdistr}}. In the univariate case, this means that \code{dp}
or \code{cp} can be vectors of longer length than indicated earlier. This
occurrence is interpreted in the sense that the additional components after
the first one are regarded as regression coefficients of a \code{selm} model,
and they are transferred unchanged to the matching components of the
transformed parameter set; the motivation is given in Section 3.1.4 of
Azzalini and Capitanio (2014). In the multivariate case, \code{dp[[1]]} and
\code{cp[[1]]} can be matrices instead of vectors; the rows beyond the first
one are transferred unchanged to \code{cp[[1]]} and \code{dp[[1]]},
respectively. 
}

\references{
Arellano-Valle, R. B. and Azzalini, A. (2013, available on-line 12 June 2011). 
The centred parameterization and related quantities of the skew-\emph{t} 
distribution. \emph{J. Multiv. Analysis} \bold{113}, 73-90. 

Azzalini, A. with the collaboration of Capitanio, A. (2014). 
\emph{The Skew-Normal and Related  Families}. 
Cambridge University Press, IMS Monographs series.

Azzalini, A. and Dalla Valle, A. (1996).
The multivariate skew-normal distribution.
\emph{Biometrika} \bold{83}, 715--726.
}

\seealso{
  \code{\link{makeSECdistr}}, \code{\link{summary.SECdistr}}, 
  \code{\link{sn.cumulants}}, 
  
  the \sQuote{Note} at \code{\link{summary.selm}} for the reason why
  \acronym{CP} is the default parameterization in that function and in 
  related ones,
  
  the \sQuote{Examples} at \code{\link{rmsn}} for use of the \acronym{CP}
  parameterization
}

\examples{
# univariate case
cp <- dp2cp(c(1, 2222, 3333, 2, 3), "SN")
dp <- cp2dp(cp, "SN")
# notice that the 2nd and the 3rd component remain unchanged
#
# multivariate case
dp3 <- list(xi=1:3, Omega=toeplitz(1/(1:3)), alpha=c(-3, 8, 5), nu=6)
cp3 <- dp2cp(dp3, "ST")
dp3.back <- cp2dp(cp3, "ST")
#
op3 <- dp2op(dp3, "ST")
dp3back <- op2dp(op3,"ST")
}

\keyword{distribution}
