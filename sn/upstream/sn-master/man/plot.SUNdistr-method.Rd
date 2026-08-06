%  file sn/man/plot.SUNdistr-method.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013-2021 Adelchi Azzalini
%---------------------
\name{plot.SUNdistr-method}
\docType{methods}
\alias{plot.SUNdistr}
\alias{plot,SUNdistr-method}
\alias{plot,SUNdistr,missing-method}
 
\title{Plotting method for class \code{SUNdistr}}

\description{Plotting method for class \code{SUNdistr}}

\usage{
\S4method{plot}{SUNdistr}(x, range, nlevels = 8, levels, npt,  main, comp, compLabs, gap = 0.5, ...)
}

\arguments{
  \item{x}{an object of class \code{SUNdistr}}

  \item{range}{in the univariate case, a vector of length 2 which defines 
    the plotting range; in the multivariate case, a matrix with two rows where 
    each column defines the plotting range of the corresponding component 
    variable. If missing, a sensible choice is made.}

  \item{nlevels}{
    number of contour levels desired \bold{iff} levels is not supplied.}
  
  \item{levels}{numeric vector of levels at which to draw contour lines.}
  
  \item{npt}{a numeric value or vector (in the univariate and in the 
    multivariate case, respectively) to assign the number of evaluation points
    of the distribution, on an equally-spaced grid over the \code{range}
    defined above.  Default value: 251 in the univariate case, a vector of
    101's in the multivariate case.}
    
  \item{main}{a character string for main title; if missing, one is built
    from the available ingredients.}
 
  \item{comp}{an optional integer vector representing the subset of the vector
    \code{1:d}, if \code{d} denotes the dimensionality of the distribution.}

  \item{compLabs}{a vector of character strings or expressions used to label
    the variables in the plot; 
    if missing, \code{slot(object,"compNames")[comp]} is used.}

  \item{gap}{a numeric value which regulates the gap between panels of a
    multivariate plot when \code{d>2}; default: \code{0.5}.}

  \item{\dots}{additional graphical parameters}
}

\details{
  For univariate density plots, \code{probs} are used to compute quantiles
  from the appropriate distribution, and these are superimposed to the plot of
  the density function, unless \code{probs=NULL}. In the multivariate case,
  each bivariate plot  is constructed as a collection of contour curves,
  one curve for each probability level; consequently, \code{probs} cannot be 
  missing or \code{NULL}. The level of the density contour lines are chosen 
  so that each curve circumscribes a region with the  quoted probability, 
  to a good degree of approssimation; for additional information, see 
  Azzalini and Capitanio (2014), specifically Complement 5.2 and p.179, 
  and references therein.}
  
\value{an invisible list. In the univariate case the list has three components:
the input object representing the distribution and two numeric vectors with 
the coordinates of the plotted density values. 
In the multivariate case, the first element of the list is the input object 
representing the distribution and all subsequent list elements are lists with 
components of the panels comprising the matrix plot; 
the elements of these sub-lists are: 
the vectors of \code{x} and \code{y} coordinates, the names of the variables, 
the density values at the \code{(x,y)} points, a vector of the density levels 
of the curves appearing in each panel plot.}

% \references{%% ~put references to the literature/web site here ~}
\author{Adelchi Azzalini} 

% \note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{makeSUNdistr}},  \code{\link{SUNdistr-class}}  }
  
 
\examples{
xi <- c(1, 0, -1)
Omega <- matrix(c(2,1,1, 1,3,1, 1,1,4), 3, 3)
Delta <- matrix(c(0.72,0.20, 0.51,0.42, 0.88, 0.94), 3, 2, byrow=TRUE)
Gamma <- matrix(c(1, 0.8, 0.8, 1), 2, 2)
dp3 <- list(xi=xi, Omega=Omega, Delta=Delta, tau=c(-0.5, 0), Gamma=Gamma)
sun3 <- makeSUNdistr(dp=dp3, name="SUN3", compNames=c("x", "w", "z"))
plot(sun3, npt=rep(51,3))
p <- plot(sun3, comp=2:3, compLabs=c(expression(x[2]), expression(x[3])))
# str(p)
}

\keyword{methods}
\keyword{hplot}
