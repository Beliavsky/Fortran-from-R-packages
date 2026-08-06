%  file sn/man/symm-modulated-distr.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2017 Adelchi Azzalini
%---------------------
\name{symm-modulated-distr}
\encoding{UTF-8}
\alias{SymmModulatedDistr}
\alias{dSymmModulated}
\alias{rSymmModulated}
\alias{dmSymmModulated}
\alias{rmSymmModulated}
\alias{plot2D.SymmModulated}
\concept{symmetry-modulated distribution}
\concept{skew-symmetric distribution}

\title{Symmetry-modulated distributions}

\description{Symmetry-modulated distributions, univariate and multivariate,
AKA skew-symmetric distributions}

\usage{
dSymmModulated(x, xi=0, omega=1, f0, G0, w, par.f0, par.G0, odd="check", 
  log=FALSE, ...)
rSymmModulated(n=1, xi=0, omega=1, f0, G0, w, par.f0, par.G0, odd="check", ...) 
dmSymmModulated(x, xi, Omega, f0, G0, w, par.f0, par.G0, odd="check", 
  log=FALSE, ...) 
rmSymmModulated(n=1, xi, Omega, f0, G0, w, par.f0, par.G0, odd="check", ...)
plot2D.SymmModulated(range, npt=rep(101,2), xi=c(0,0), Omega, f0, G0, w, 
  par.f0, par.G0, odd="check", ...)  
}

\arguments{
  \item{x}{a vector of coordinates where the density must be evaluated; 
   for multivariate densities, evaluated by \code{dmSymmModulated}, 
   a matrix is also allowed, each row representing a point.}
   
  \item{xi}{a numeric vector representing the location parameter; 
   if must have length 1 for \code{dSymmModulated} and \code{rSymmModulated},
   length 2 for \code{plot2D.SymmModulated}).}
   
  \item{omega}{a positive value representing the scale parameter.}
  
  \item{f0}{a character string denoting the symmetric density to be modulated;
    admissible values for \code{dSymmModulated} and \code{dSymmModulated} 
    are \code{"beta"}, \code{"cauchy"}, \code{"logistic"} (or \code{"logis"}), 
    \code{"normal"} (or \code{"norm"}), \code{"t"}, \code{"uniform"};   
    for the other functions the possible 
    values are \code{"cauchy"}, \code{"normal"} (or \code{"norm"}), \code{"t"}; 
    the meaning of the names is described in the \sQuote{Details}  section.}
   
  \item{G0}{a character string denoting the symmetric distribution used in the
    modulating factor; admissible values are \code{"beta"}, \code{"cauchy"}, 
    \code{"logistic"} (or \code{"logis"}), \code{"normal"} (or \code{"norm"}), 
    \code{"t"}, \code{"uniform"},  with meaning  described in the
    \sQuote{Details} section.}
    
  \item{w}{the name (\emph{not} as a character string) of a user-defined 
    function which satisfies the condition \eqn{w(-z)=-w(z)} for all \eqn{z};
    see the \sQuote{Details} section for additional specifications.}    
 
  \item{par.f0, par.G0}{parameters required by \code{f0} and \code{G0}, 
    when they are of type \code{"beta"} or \code{"t"}, otherwise ignored.}
   
  \item{odd}{a character string, with possible values \code{"check"} (default),
    "assume", "force",  for regulation of the behaviour about the condition
    that \code{w} is an odd function, as explained in the \sQuote{Details}
    section.} 
 
  \item{log}{logical (default: \code{FALSE}); 
    if \code{TRUE}, densities are given as log(densities).} 

  \item{n}{an integer value (default: \code{n=1}) indicating the number of 
    random numbers.}
   
  \item{Omega}{a symmetric positive-definite matrix which regulates the 
    dependence structure of \code{f0} and so of the final density.}
  
  \item{range}{a two-column matrix whose column-wise range is taken as the
    plotting intervals  on the coordinated axes forming a bivariate grid of
    points over which the density is plotted.}

  \item{npt}{a numeric vector with two elements representing the number of  
    equally-spaced points on each axis spanning the \code{range} described 
    above;  default value is \code{rep(101,2)}.}
    
  \item{\dots}{optional parameters regulating the function \code{w} and, for 
    \code{plot2D.SymmModulated} only, graphical parameters to be supplied to
    function \code{\link[graphics]{contour}}. }
  
}

\section{Background}{
In the univariate case, start from symmetric density function \eqn{f_0}{f₀}, 
such that  \eqn{f_0(z)=f_0(-z)}{f₀(z)=f₀(-z)} for all \eqn{z}, and 
\sQuote{modulate}  it in the form 
\deqn{f(z) = 2\, f_0(z)\, G_0\{w(z)\}}{f(z) = 2 f₀(z) G₀\{w(z)\} }
where \eqn{G_0}{G₀} is a univariate symmetric (about 0) distribution function 
and \eqn{w(z)}is a real-valued odd function, hence satisfying the condition 
\eqn{w(-z)=-w(z)};
then $f(z)$ is a proper density function wich integrates to 1. 
A subsequent location and scale transformation applied to \eqn{f(z)} 
delivers the final density. 
Specifically, if \eqn{Z} denotes a univariate random variable with density 
\eqn{f(z)}, then the computed density pertains to the transformed variable
\deqn{\xi + \omega Z.}

In the multivariate case, the scheme is similar, with natural adaptation.
Density \eqn{f_0}{f₀} is  now \eqn{d}-dimensional, while \eqn{G_0}{G₀} 
is still univariate. The conditions \eqn{f_0(z)=f_0(-z)}{f₀(z)=f₀(-z)} 
and \eqn{w(-z)=-w(z)} refer to a \eqn{d}-dimensional vector \eqn{z}.  
Given a \eqn{d \times d}{d x d} symmetric positive-definite matrix 
\eqn{\Omega}, we extract the the square roots \eqn{\omega} of the diagonal 
element of  \eqn{\Omega} and correspondingly obtain the scale-free matrix
\deqn{\bar\Omega = \mathrm{diag}(\omega)^{-1}\, \Omega\,
                             \mathrm{diag}(\omega)^{-1}
     }{cov(\Omega) = diag(\omega)⁻¹ \Omega diag(\omega)⁻¹}
which is used to regulate the dependence structure of \eqn{f_0(z)}{f₀(z)} 
and so of \eqn{f(z)}. 
If \eqn{Z} is multivariate random variable with density \eqn{f(z)}, 
then the final distribution refers to 
\deqn{\xi + \mathrm{diag}(\omega)\,Z }{\xi + diag(\omega) Z }
where \eqn{\xi} is a \eqn{d}-dimensional vector of location parametes.

This construction was put forward by Azzalini and Capitanio (2003). 
An essentially equivalent formulation has been presented by Wang et al. (2004).
A summary account is available in Section 1.2 of Azzalini and Capitanio
(2014); this includes, inter alia, an explanation  of why the term
\sQuote{symmetry-modulated} distributions is preferred to 
\sQuote{skew-symmetric} distributions. 

Random number generation is based on expression (1.11a) of
Azzalini and Capitanio (2014).
}

\section{Details}{
Functions \code{dSymmModulated} and \code{rSymmModulated} deal with univariate 
distributions, for computing densities and generating random rumbers, 
respectively; 
\code{dmSymmModulated} and \code{rmSymmModulated} act similarly 
for multivariate distributions. For the bivariate case only,
\code{plot2D.SymmModulated} computes a density over a grid of
coordinates and produces a \code{contour} plot. 
 
The distribution names used in \code{f0} and \code{G0} have, in the univariate
case, the same meaning as described in the \code{\link[stats]{Distributions}}
page, with the following exceptions, to achive symmetry about 0: 
\code{"uniform"} denotes a  uniform distribution over the interval 
\eqn{(-1, 1)}; \code{"beta"} denotes the a symmetric Beta distribution with
support  over the interval \eqn{(-1, 1)} and a common value of the shape 
parameters.

In the multivariate case, the available options \code{"normal"} and \code{"t"}
for \code{f0} refer to densities computed by  \code{\link[mnormt]{dmnorm}} 
and \code{\link[mnormt]{dmt}} with 0 location and correlation matrix   
\eqn{\bar\Omega}{cor(\Omega)}, implied by \eqn{\Omega}. 
Argument \code{G0} has the same meaning as in the univariate case.

Options \code{"beta"} and \code{"t"} for \code{f0} and \code{G0} require the
specification of a shape parameter, via the arguments \code{par.f0} and
\code{par.G0}, respectively. For \code{"beta"} the parameter represents
the common value of the shape parameters of \code{\link[stats]{Beta}};
for \code{"t"}, it represents \code{df} of \code{\link[stats]{TDist}} and
\code{\link[mnormt]{dmt}}. 
 
Function \code{w} most be  of the form \code{w <- function(z, ...)} where
\code{...} are optional additional parameters and \code{z} represents valued
of the standardized form of the density; in the univariate case, \code{x} and
\code{z} are related by \code{z=(x-xi)/omega} and an analogous fact holds in
the multivariate setting. The function must satisfy the condition
\eqn{w(-z)=-w(z)}. It is assumed that the function is vectorized and, in the
multivariate case, it will be called with \code{z} representing a matrix with
\code{d} columns, if \code{d} denotes the dimensionality of the random
variable.
 
Argument \code{odd} regulates the behaviour with respect to the condition 
\eqn{w(-z)=-w(z)}. If its value is \code{"assume"}, the condition is just
assumed to hold, and no action is taken. If the value is \code{"check"}
(deafult), a \emph{limited} check is performed; namely, in case of densities,
the check is at 0 and the supplied \code{x} points, while for random numbers  
the check is at 0 and the generated points.
The value \code{"force"} ensures that the condition is satisfied by 
actually constructing a modified version of the user-supplied function 
\code{w}, such that the required condition is enforced. 
}

\value{
For \code{dSymmModulated}, \code{rSymmModulated} and  \code{dmSymmModulated}, 
a numeric vector; for \code{dmSymmModulated} a matrix, unless \code{n=1}.

For \code{plot2D.SymmModulated} an invisible list containing the \code{x}
and \code{y} coordinates forming the grid over which the density \code{pdf}
has been evaluated for plotting.
} 

\author{Adelchi Azzalini} 
 
\references{

Arellano-Valle, R. B., Gómez, H. W. and Quintana, F. A. (2004).
 A new class of skew-normal distributions. \emph{Comm. Stat., Theory &
 Methods}, \bold{58}, 111-121. 

Azzalini, A. and Capitanio, A. (2003).
 Distributions generated by perturbation of symmetry 
 with emphasis on a multivariate skew-\emph{t} distribution.
 \emph{J.Roy. Statist. Soc. B} \bold{65}, 367--389.
 Full version of the paper at \url{https://arXiv.org/abs/0911.2342}
 
Azzalini, A. with the collaboration of Capitanio, A. (2014). 
 \emph{The Skew-Normal and Related Families}. 
 Cambridge University Press, IMS Monographs series.    

Wang, J., Boyer, J. and Genton, M. G. (2004).
 A skew-symmetric representation of multivariate distributions.
 \emph{Statistica Sinica}, \bold{14}, 1259-1270.
}

\seealso{
\code{\link[stats]{Distributions}}, 
\code{\link[stats]{Beta}}, \code{\link[stats]{TDist}},
\code{\link[mnormt]{dmnorm}}, \code{\link[mnormt]{dmt}}, 
\code{\link[graphics]{contour}}
}


\examples{
x <- seq(2, 13, length=45)
wLinear <- function(z, lambda) lambda*z
y <- dSymmModulated(x, 5, 2, f0="normal", G0="normal", w=wLinear, lambda=3)
# the same of dsn(x, 5, 2, 3), up to negligible numerical differences
#
wSGN <- function(z, lambda) z*lambda[1]/sqrt(1 + lambda[2]*z^2) 
y <- dSymmModulated(x, 5, 2, f0="normal", G0="normal", w=wSGN, lambda=c(3,5))
# SGN distribution of Arellano-Valle et al. (2004)
#
wST <- function(z, lambda, nu) lambda*z*sqrt((nu+1)/(nu+z^2))
y <- rSymmModulated(n=100, 5, 2, f0="t", G0="t", w=wST, par.f0=8, par.G0=9, 
        lambda=3, nu=8)
# equivalent to rst(n=100, 5, 2, 3, 8) 
#
wTrigs <- function(z, p, q) sin(z * p)/(1 + cos(z * q))
x <- seq(-1, 1, length=51)
y <- dSymmModulated(x, 0, 1, f0="beta", G0="logistic", w=wTrigs, par.f0=2,
        par.G0=NULL, p=5, q=0.5)   
plot(x, y, type="l")        
# univariate analogue of the bivariate distribution on pp.372-3 of
# Azzalini & Capitanio (2003)           
#
range <- cbind(c(-3,3), c(-3,3))
wMvTrigs <- function(z, p, q) sin(z \%*\% p)/(1 + cos(z \%*\% q))
plot2D.SymmModulated(range, xi=c(0,0), Omega=diag(2), f0="normal", G0="normal",
    w=wMvTrigs, par.f0=NULL, par.G0=NULL, p=c(2,3), q=c(1,1), col=4)
# w(.) as in (1.6) of Azzalini & Capitanio (2014, p.4) and plot as in 
# bottom-right panel of their Figure 1.1.      
}
    
\keyword{distributions}
\keyword{hplot}
