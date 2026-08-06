%  file sn/man/pprodt2.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013-2020 Adelchi Azzalini
\name{pprodt2}
\alias{pprodt2}
\alias{qprodt2}
\alias{pprodn2}
\concept{bivariate normal distribution}
\concept{bivariate Student's t distribution}
\concept{product of two random variables} 

\title{
The distribution function of the product of two jointly normal or Student's \emph{t} variables
}
\description{
Consider the product \eqn{W=X_1 X_2} from a bivariate random variable 
\eqn{(X_1, X_2)} having joint normal or Student's \emph{t} distribution, 
with 0 location and unit scale parameters. 
Functions are provided for the distribution function of \eqn{W} in the 
normal and the \emph{t} case, and the quantile function for the \emph{t} case.
}
\usage{
pprodn2(x, rho)
pprodt2(x, rho, nu)
qprodt2(p, rho, nu, tol=1e-5, trace=0)
}
\arguments{
  \item{x}{a numeric vector}
  \item{p}{a numeric vector of probabilities}
  \item{rho}{a scalar value representing the correlation (or the matching
         term in the \emph{t} case when correlation does not exists)}
  \item{nu}{a positive scalar representing the degrees of freedom}
  \item{tol}{the desired accuracy (convergence tolerance),
        passed to function \code{\link[stats]{uniroot}} }
  \item{trace}{integer number for controlling tracing information, 
       passed on to \code{uniroot}}
}
\details{
Function \code{pprodn2} implements results for the central case in
Theorem 1 of Aroian \emph{et al.} (1978).

Function \code{pprodt2} implements formulae in Theorem 1 of Wallgren (1980).
Corresponding quantiles are obtained by \code{qprodt2} by solving the 
pertaining non-linear equations with the aid of \code{\link[stats]{uniroot}},
one such equation for each element of \code{p}. 
}
\value{a numeric vector}
\references{
Aroian, L.A., Taneja, V.S,  &  Cornwell, L.W. (1978).
Mathematical forms of the distribution of the product of two normal variables.
\emph{Communications in statistics. Theory and methods}, 7, 165-172

Wallgren, C. M. (1980).
The distribution of the product of two correlated \emph{t} variates.
\emph{Journal of the American Statistical Association}, 75, 996-1000
}
\author{Adelchi Azzalini}
% \note{ ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link[stats]{uniroot}}  }
\examples{
p <-  pprodt2(-3:3, 0.5, 8)
qprodt2(p, 0.5, 8)
}
%
\keyword{distribution} % use one of  RShowDoc("KEYWORDS")

