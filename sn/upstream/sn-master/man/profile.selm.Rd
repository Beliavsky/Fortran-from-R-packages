%  file sn/man/profile.selm.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2015-2026 Adelchi Azzalini
%---------------------

\name{profile.selm}
\alias{profile.selm}
\alias{profile,selm-method}      
\concept{profile log-likelihood function}
\concept{confidence interval}
\concept{confidence region}

\title{Profile log-likelihood function of selm-class objects}

\description{
One- or two-dimensional profile (penalized) log-likelihood function 
of a \code{selm} fit and corresponding confidence interval or regions.
}
\usage{
  \S4method{profile}{selm}(fitted, param.type, param.name, param.values, npt, 
   opt.control = list(), plot.it = TRUE, log = TRUE, levels, trace = FALSE, ...)
}

\arguments{
  \item{fitted}{an object of class \code{selm} as produced by a call to
   function \code{selm} with a univariate response.}
  
  \item{param.type}{
   a character string with the required parameterization; it must be either
   \code{"CP"} or \code{"DP"}, or possibly their equivalent lowercase.}
   
  \item{param.name}{
   either a single character string  or a vector of two such terms with the 
   name(s) of the parameter(s) for which the profile log-likelihood is 
   required; these names must match those appearing in 
   \code{\link{summary.selm}(object, param.type)}.}
  
  \item{param.values}{in the one-parameter case, a numeric vector with 
   the values where the log-likelihood must be evaluated; in the two-parameter
   case, a list of two such vectors used to build a  grid of coordinates 
   of points. Their range must identify an interval or a rectangle which
   includes the \acronym{MLE} or \acronym{MPLE} obtained by \code{\link{selm}}.
   See \sQuote{Details} for more information.}
  
  \item{npt}{
   in case the vector or any of the vectors of argument \code{param.values} 
   has length 2, an equally spaced grid of values is build with length equal
   to the corresponding component of \code{npt}.
   If the above condition is met but this argument is missing,  
   a default choice is made, namely 51 or (26,26) in the one- or two-parameter
   case, respectively.}
  
  \item{opt.control}{
   an optional list passed as argument \code{control} to \code{optim} to
   optimize the log-likelihood; see \sQuote{Details} for more information.}
   
  \item{plot.it}{a logical value; if \code{TRUE} (default value), a plot 
   is produced representing the deviance, which is described in 
   \sQuote{Details} below.
   In the one-parameter case, a confidence interval of prescribed \code{level}
   is marked on the plot; in the two-parameter case, the contour curves are
   labelled with approximate confidence levels. See however 
   for more information.
   }
  
  \item{log}{a logical value (default: \code{TRUE})  indicating whether the
   scale and tail-weight parameter (the latter only for the \acronym{ST} 
   family) must be log-transformed, if case any of them occurs in   
   \code{param.name}. This applies to \code{omega} and 
   \code{nu} in the \code{DP} parameter set and to \code{s.d.} and 
   \code{gamma2} in the \code{CP} parameter set.}
   
  \item{levels}{a single probability value (in the one-parameter case) or a 
   vector of such values (in the two-parameter case) for which the confidence
   interval or region is requited; if missing, the default value is \code{0.95}
   in the one-parameter case and \code{c(0.25, 0.5, 0.75, 0.9, 0.95, 0.99)}
   in the two-parameter case.}
   
  \item{trace}{a logical value (default: \code{FALSE}) to activate printing
   of intermediate outcome of the log-likelihood optimization process}
  
  \item{\dots}{optional graphical parameters passed to the plotting functions.}  
}

\details{
For each chosen point of the parameter(s) to be profiled, the 
log-likelihood is maximized with respect to the remaining parameters.
The optimization process is accomplished using the \code{\link[stats]{optim}}
optimization function, with \code{method="BFGS"}. This step can be regulated by
the user via \code{opt.control} which is  passed to \code{\link[stats]{optim}} 
as \code{control} argument, apart from element \code{fnscale} whose use is
reserved.

If the original \code{fitted} object included a fixed parameter value, 
this is kept fixed here. If the estimation method was \code{"MPLE"},
that choice carries on here; in case the penalty function was user-defined,
it must still be accessible.

For plotting purposes and also in the numerical output, the deviance
function \eqn{D} is used, namely
\deqn{D = 2\left[\max(\log L) - \log L\right]}{D = 2*[max(log L)- log L]}
where \eqn{L} denotes the likelihood.  

The range of \code{param.values} must enclose the maximum (penalized)
likelihood estimates (\acronym{MLE} or \acronym{MPLE}) by an adequate extent 
such that suitable confidence intervals or regions can be established from
standard asymptotic theory. If this condition does not hold, the function
still proceeds, but no confidence interval or region is delivered.
Appropriate ranges of the parameters can be selected with the help of a 
preliminary call to \code{\link[=confint.selm]{confint}}.

For the \acronym{SN} family and \acronym{DP} parameterization, the asymptotic
theory is actually non-standard near the important point \eqn{\alpha=0}, 
but the correspondence with the regular case of the \acronym{CP} 
parameterization, still allows to derive confidence regions using standard
procedures; for additional information, see Section 3.1.6 of
Azzalini and Capitanio (2014).
When the \acronym{MLE} occurs on the frontier of the parameter space, 
a message is issued and no confidence interval is produced, while in the
two-parameter case the plot is not labelled with probability values, but only
with deviance levels.
}

\value{An invisible list whose components, described below, are partly 
different in the one- and the two-parameter cases.
  \item{call}{the calling statement}
  \item{<param1>}{values of the first parameter}
  \item{<param2>}{values of the second parameter (in the two-parameter case)}
  \item{logLik}{numeric vector or matrix of the profile log-likelihood values}
  \item{confint}{in the one-parameter case, the confidence interval}
  \item{level}{in the one-parameter case, the confidence level}
  \item{deviance.contour}{in the two-parameter case, a list of lists whose
    elements identify each curve of the contour plot}
  }
   

\author{Adelchi Azzalini}

\references{
 Azzalini, A. with the collaboration of Capitanio, A. (2014). 
 \emph{The Skew-Normal and Related Families}. 
 Cambridge University Press, IMS Monographs series. 
}
 
% \note{}

\section{Warnings}{
\itemize{
\item This function is experimental and changes in future versions
of the package may occur. Users should not rely on the persistence of the 
same user interface or the same name(s).

\item It is a known fact that, in some critical situations, peculiar outcomes
are produced.
}} 

\seealso{
\code{\link{selm}} and related \sQuote{methods} for the 
  \code{\linkS4class{selm}-class}, notably \code{\link[=confint.selm]{confint}}.

\code{\link{makeSECdistr}} for the \acronym{CP}/\acronym{DP} parameterizations,

\code{\link[stats]{optim}} for its \code{control} argument
}

\examples{
data(ais, package="sn")
m1 <-  selm(log(Fe) ~ BMI + LBM, family = "sn", data = ais)

pll <- profile(m1, "dp", param.name="alpha", param.val=c(-3,2))

profile(m1, "cp", param.name="gamma1", param.val=seq(-0.7, 0.4, by=0.1))

# in the next example, the grid points are reduced to save execution time
pll <- profile(m1, "cp", param.name=c("(Intercept.CP)", "gamma1"),
           param.val = list(c(1.5, 4), c(-0.8, 0.5)), npt=c(11,11) )
}

