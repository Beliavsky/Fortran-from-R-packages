%  file sn/man/SUNdistr-vase.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013-2023 Adelchi Azzalini
%---------------------
\name{SUNdistr-base}
\encoding{UTF-8}

\alias{SUNdistr-base}
\alias{dsun}
\alias{psun}
\alias{rsun}
\alias{sunMean}
\alias{sunVcov}
\alias{sunMardia}

\title{The Unified Skew-Normal (SUN) probability distribution}

\description{
Density, distribution function, random number generation, the mean value,
the variance-covariance matrix and the Mardia's measures  of multivariate
skewness and kurtosis of the \acronym{SUN} probability distribution.
}

\usage{
dsun(x, xi, Omega, Delta, tau, Gamma, dp = NULL, log = FALSE, silent=FALSE, ...)
psun(x, xi, Omega, Delta, tau, Gamma, dp = NULL, log = FALSE, silent=FALSE, ...)
rsun(n=1, xi, Omega, Delta, tau, Gamma, dp = NULL, silent=FALSE)
sunMean(xi, Omega, Delta, tau, Gamma, dp = NULL, silent=FALSE, ...)
sunVcov(xi, Omega, Delta, tau, Gamma, dp = NULL, silent=FALSE, ...)
sunMardia(xi, Omega, Delta, tau, Gamma, dp = NULL, silent=FALSE, ...)
}
\arguments{
  \item{x}{either a vector of length \code{d}, where \code{d=ncol(Omega)},
   with the coordinates of the point where the density or the
   distribution function must be evaluated, or alternatively
   a \code{d}-column matrix whose rows represent a set of points.}

  \item{xi}{a numeric vector of length \code{d} representing the 
   location parameter of the distribution; see \sQuote{Background}.
   In a call to \code{dsun} and \code{psun}, \code{xi} can be a matrix,
   whose rows represent a set of location parameters;
   in this case, its dimensions must match those of \code{x}.}

  \item{Omega}{a symmetric positive definite matrix of dimension \code{(d,d)};
   see \sQuote{Details}.}

  \item{Delta}{a matrix of size \code{(d,m)}, where \code{m=length(tau)};
   see \sQuote{Details} about its constraints.}

  \item{tau}{a vector of length \code{m}, say.}
  
  \item{Gamma}{a symmetric positive definite matrix of dimension \code{(m,m)}
   with 1's on its main diagonal, that is, a correlation matrix}
    
  \item{dp}{a list with five elements, representing \code{xi} 
   (which must be a vector in this case), 
   \code{Omega}, \code{Delta}, \code{tau} and \code{Gamma}, 
   with restrictions indicated in the \sQuote{Details}.
   Its default value is \code{NULL};  if \code{dp} is assigned, 
   the individual parameters must not be specified.}
   
  \item{n}{a positive integer value.} 
    
  \item{log}{a logical value (default value: \code{FALSE}); 
   if \code{TRUE}, log-densities and log-probabilities are returned.}
  
  \item{silent}{a logical value which indicates the action to take in the case
   \code{m=1}, which could be more convenently handled by functions for the
   \acronym{SN/ESN} family. If \code{silent=FALSE} (default value), a warning
   message is issued; otherwise this is suppressed.}
    
  \item{\dots}{additional tuning arguments passed either to 
    \code{\link[mnormt]{pmnorm}} (for \code{dsun}, 
    \code{psun} and \code{sunMean}) or to \code{\link[mnormt]{mom.mtruncnorm}} 
    (for \code{sunVcov} and \code{sunMardia}); see also \sQuote{Details}. 
    }
}

\details{
A member of the \acronym{SUN} family of distributions is identified by five 
parameters, which are described  in the \sQuote{Background} section.  
The five parameters can be supplied by combining them in a list, denoted
\code{dp}, in which case the individual parameters must \emph{not} be
supplied. The elements of \code{dp} must appear in the above-indicated order
and must be named.

The optional arguments in \code{...} are passed to \code{\link[mnormt]{pmnorm}},
which uses \code{\link[mnormt]{ptriv.nt}} when \code{d=3}, 
\code{\link[mnormt]{biv.nt.prob}} when \code{d=2} and 
and \code{\link[mnormt]{sadmvn}} when \code{d>2}.
In practice these arguments are effective only if \code{d>3}, 
since for lower dimensions the computations are made to full available 
precision anyway.
A similar fact applies to the \code{...} argument passed to 
\code{\link[mnormt]{mom.mtruncnorm}}. 

Some numerical inaccuracy is inevitably involved in these computations.
In most cases,  they are of negligible extent, but they can possibly become
more relevant, especially in the computation of higher order moments
involved by \code{sunMardia}, depending on the dimension \code{d} and on 
the specific parameter values.
Consider the \sQuote{Warning} section in \code{\link[mnormt]{recintab}}
which is used by \code{\link[mnormt]{mom.mtruncnorm}}. 

Note that the Mardia measures are computed in their \eqn{\gamma}-form, not the 
\eqn{\beta}-form; in practice, these only differ by a \eqn{d(d+2)} shift 
in the measure of kurtosis, if \eqn{d} denotes the dimension.

The above-described functions operate following the traditional \R scheme for 
probability distributions. Another scheme, coexisting with the classical one, 
works with \code{SUNdistr-class} objects, which represent \acronym{SUN} 
distributions, by encapsulating their parameters and other characteristics. 
These objects are created by \code{\link{makeSUNdistr}}, and various methods 
exist for them; see \code{\link{SUNdistr-class}}. 
Moreover these objects can be manipulated by a number of tools, described 
in \code{\link{SUNdistr-op}}, leading to new objects of the same class.

Given the constraints on the terms \code{Omega, Delta, Gamma}, 
\emph{it is advisable}
to always use the function \code{\link{makeSUNdistr}} to check the joint 
admissibility of the component parameters, \emph{even if} one is not going
to actually use the returned object and the related methods and functions. 
Consider that, for efficiency reasons, the functions following the 
traditional \R scheme for probability distributions  
(\code{dsun}, ..., \code{sunMean}, ...)
do not perform those admissibility checks.
}

\value{
The structure of the returned value depends on the called function, as follows:
 \tabular{rl}{
 \code{dsun, psun}  \tab a vector of length \code{nrow(x)} representing 
                          density or probability values, \cr
             \tab or their log-transformed values if \code{log=TRUE},\cr
 \code{rsun} \tab a matrix of size \code{(n,d)}, 
                  where each row represents a \acronym{SUN} random vectors,\cr
 \code{sunMean} \tab a vector of length \code{d} representing the mean value,\cr
 \code{sunVcov} \tab a matrix of size \code{(d,d)}  representing the 
   variance-covariance matrix,\cr
 \code{sunMardia} \tab a vector of length two with the Mardia's measures of
   multivariate skewness and kurtosis.
 }
}

\section{Background}{
A member of the \acronym{SUN} family is characterized by two
dimensionality indices, denoted \eqn{d} and \eqn{m}, and a set of five
parameters blocks (vector and matrices, as explained soon).  
The value \eqn{d} represents the number of observable
components; the value \eqn{m} represents the number of latent (or hidden)
variables notionally involved in the construction of the distribution.
The parameters and their corresponding \R variables are as follows:
 \tabular{rcl}{
 \eqn{\xi}    \tab \code{xi}    \tab a vector of length \eqn{d}, \cr
 \eqn{\Omega} \tab \code{Omega} \tab a matrix of size \eqn{(d,d)}, \cr
 \eqn{\Delta} \tab \code{Delta} \tab a matrix of size \eqn{(d,m)}, \cr
 \eqn{\tau}   \tab \code{tau}   \tab a vector of length \eqn{m}, \cr
 \eqn{\Gamma} \tab \code{Gamma} \tab a matrix of size \eqn{(m,m)},  
 }
and must satisfy the following  conditions:  
\enumerate{
\item \eqn{\Omega} is a symmetric positive definite matrix;
\item \eqn{\Gamma} is a symmetric positive definite matrix with 1's 
    on the main diagonal, hence a correlation matrix;
\item if \eqn{\bar\Omega}  % {\Omega°} 
  denotes the correlation matrix associated to \eqn{\Omega},  
  the matrix of size \eqn{(d+m)\times(d+m)}{((d+m), (d+m))}
  \if{html}{formed by the \eqn{2 \times 2} blocks}
  \ifelse{latex}{
         \deqn{\left(\begin{array}{cc} 
            \bar\Omega  & \Delta \\
            \Delta^\top & \Gamma
            \end{array} \right)} 
       }{ % non-LaTeX
        \tabular{rrcc}{
            \tab  \tab \eqn{\bar\Omega}  \tab  \eqn{\Delta} \cr
            \tab  \tab \eqn{\Delta'}  \tab  \eqn{\Gamma}  }
       }
  must be a positive definite correlation matrix.}
 
The formulation adopted here has arisen as the evolution of earlier 
constructions, which are recalled very briefly next.
A number of extensions of the multivariate skew-normal distributions, 
all involving a number \code{m} (with \eqn{m\ge1}) of latent variables 
(instead of \code{m=1} like the skew-normal distribution), 
have been put-forward in close succession in the years 2003-2005. 
Special attention has been drawn by the \sQuote{closed skew-normal (CSN)}
distribution  developed by González-Farías \emph{et alii} (2004a, 2004b)
and the \sQuote{fundamental skew-normal  (FUSN)} distribution
developed by Arellano-Valle and Genton (2005),
but other formulations have been considered too.

Arellano Valle and Azzalini (2006) have shown the essential equivalence
of these apparently alternative constructions, after appropriate
reparameterizations, and underlined the necessity of removing 
over-parameterizations in some cases, to avoid lack of identifiability. 
This elaboration has led to the \acronym{SUN} formulation. 
A relatively less technical account of their development is provided 
in Section 7.1 of Azzalini and Capitanio (2014), using very slightly 
modified notation and parameterization, which are the ones adopted here.  

Additional results have been presented by Arellano-Valle 
and Azzalini (2021), such as expressions for the variance matrix
and higher order moments, the Mardia's measures of multivariate  skewness
and kurtosis, which are implemented here. Another result is the
conditional distribution when the conditioning event is represented
by an orthant.
}


\references{
Arellano-Valle, R. B., and Azzalini, A. (2006).
On the unification of families of skew-normal distributions.
\emph{Scand. J. Stat.} \bold{33}, 561-574. 
Corrigendum in  \bold{49} (2022), 1418-1419.

Arellano-Valle, R. B. and Azzalini, A. (2021).
Some properties of the unified skew-normal distribution.
\emph{Statistical Papers} \bold{63}, 461-487,
\doi{https://doi.org/10.1007/s00362-021-01235-2}.
Corrigendum in \bold{65} (2024), 1133, 
\doi{https://doi.org/10.1007/s00362-023-01412-5}
% see also \href{https://arxiv.org/abs/2011.06316}{arXiv:2011.06316}
  
Arellano-Valle, R. B. and Genton, M. G. (2005). 
On fundamental skew distributions. 
\emph{J. Multivariate Anal.} \bold{96}, 93–1116.

Azzalini, A. with the collaboration of Capitanio, A. (2014). 
\emph{The Skew-Normal and Related Families}. 
Cambridge University Press, IMS Monographs series.

González-Farías, G., Domínguez-Molina, J. A., & Gupta, A. K. (2004a). 
Additive properties of skew normal random vectors.
\emph{J. Statist. Plann. Inference} \bold{126}, 521-534.

González-Farías, G., Domínguez-Molina, J. A., & Gupta, A. K. (2004b). 
The closed skew-normal distribution. 
In M. G. Genton (Ed.), \emph{Skew-elliptical Distributions and Their 
Applications: a Journey Beyond Normality}, Chapter 2, (pp. 25–42). 
Chapman & Hall/\acronym{CRC}. 
}

\author{Adelchi Azzalini}

\note{
The present structure and user interface of this function, and of other ones 
related to the \acronym{SUN} distribution, must be considered experimental, 
and they might possibly change in the future.}

\seealso{
  \code{\link{makeSUNdistr}} to build a \acronym{SUN} distribution object, 
  with related methods in \code{\link{SUNdistr-class}},
  and other facilities in \code{\link{SUNdistr-op}}
  
  \code{\link{convertCSN2SUNpar}} to convert a parameter set of the Closed 
  Skew-Normal formulation to the equivalent \acronym{SUN} parameter set
}

\examples{
xi <- c(1, 0, -1)
Omega <- matrix(c(2,1,1, 1,3,1, 1,1,4), 3, 3)
Delta <- matrix(c(0.72,0.20, 0.51,0.42, 0.88, 0.94), 3, 2, byrow=TRUE)
Gamma <- matrix(c(1, 0.8, 0.8, 1), 2, 2)
dp3 <- list(xi=xi, Omega=Omega, Delta=Delta, tau=c(-0.5, 0), Gamma=Gamma)
x <- c(0.8, 0.5, -1.1)
f1 <- dsun(x, xi, Omega, Delta, c(-0.5, 0), Gamma) # mode 1
f2 <- dsun(x, dp=dp3)   # mode 2, equivalent to mode 1
set.seed(1)
xm <- rsun(10, dp=dp3)
f3 <- dsun(xm, dp=dp3) 
psun(xm, dp=dp3)
sunMean(dp=dp3)
sunVcov(dp=dp3)
sunMardia(dp=dp3)
}

\keyword{distribution}
\keyword{multivariate} 
\concept{SUN distribution}
\concept{Unified Skew-Normal distribution}
\concept{CSN distribution}
\concept{Closed Skew-Normal distribution}
\concept{FUSN distribution}
\concept{Fundamental Skew-Normal distribution}
