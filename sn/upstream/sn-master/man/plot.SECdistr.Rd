%  file sn/man/plot.SECdistr.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013 Adelchi Azzalini
%---------------------
\name{plot.SECdistr}
\docType{methods}
%\alias{plot,ANY,ANY-method}
% \alias{plot,profile.mle,missing-method}
% \alias{show,SECdistrMv-method}
% \alias{show,SECdistrUv-method}
\alias{plot.SECdistr}
\alias{plot.SECdistrUv}
\alias{plot.SECdistrMv}
\alias{plot,SECdistrMv,missing-method}
\alias{plot,SECdistrUv,missing-method}
\alias{plot,SECdistrMv-method}
\alias{plot,SECdistrUv-method}

\title{Plotting methods for classes \code{SECdistrUv} and \code{SECdistrMv}}

\description{Plotting methods for classes \code{SECdistrUv} and \code{SECdistrMv}}

\usage{
\S4method{plot}{SECdistrUv}(x, range, probs, main, npt = 251, \dots)

\S4method{plot}{SECdistrMv}(x, range, probs, npt, landmarks = "auto",
   main, comp, compLabs, data = NULL, data.par = NULL, gap = 0.5, \dots) 
}

\arguments{
  \item{x}{an object of class \code{SECdistrUv} or \code{SECdistrMv}.}
   
  % \item{y}{not used, required by the generic \code{plot(x, y, ...)} function.}
   
  \item{range}{in the univariate case, a vector of length 2 which defines 
    the plotting range; in the multivariate case, a matrix with two rows where 
    each column defines the plotting range of the corresponding component 
    variable. If missing, a sensible choice is made.}

  \item{probs}{a vector of probability values. In the univariate case, the
    corresponding quantiles are plotted on the horizontal axis; it can be
    skipped by setting \code{probs=NULL}. In the multivariate case, each 
    probability value corresponds to a contour level in each bivariate plot; 
    at least one probability value is required. See \sQuote{Details} for
    further information. Default value: \code{c(0.05, 0.25, 0.5, 0.75, 0.95)}
    in the univariate case, \code{c(0.25, 0.5, 0.75, 0.95)} in the
    multivariate case.}
    
  \item{npt}{a numeric value or vector (in the univariate and in the 
    multivariate case, respectively) to assign the number of evaluation points
    of the distribution, on an equally-spaced grid over the \code{range}
    defined above.  Default value: 251 in the univariate case, a vector of
    101's in the multivariate case.}
    
  \item{landmarks}{a character string which affects the placement of some
    landmark values in the multivariate case, that is, the origin, the mode
    and the mean (or its substitute pseudo-mean), which are all aligned.
    Possible values: \code{"proper"}, \code{"pseudo"}, \code{"auto"}
    (default), \code{""}. The option \code{""} prevents plotting of the
    landmarks. With the other options, the landmarks are plotted, with some
    variation in the last one: \code{"proper"} plots the proper mean value,
    \code{"pseudo"} plots the pseudo-mean, useful when the proper mean does
    not exists, \code{"auto"} plots the proper mean if it exists, otherwise it
    switches automatically to the pseudo-mean. See \code{\link{dp2cp}} for
    more information on pseudo-\acronym{CP} parameters, including pseudo-mean.}
    
  \item{main}{a character string for main title; if missing, one is built
    from the available ingredients.}
  
  \item{comp}{a subset of the vector \code{1:d}, if \code{d} denotes the
    dimensionality of the multivariate distribution.}
    
  \item{compLabs}{a vector of character strings or expressions used to denote
    the variables in the plot; 
    if missing, \code{slot(object,"compNames")} is used.}
    
  \item{data}{an optional set of data of matching dimensionity of
    \code{object} to be superimposed to the plot. 
    The default value \code{data=NULL} produces no effect. 
    In the univariate case, data are plotted using \code{\link[graphics]{rug}} 
    at the top horizontal axis, unless if \code{probs=NULL}, in which case 
    plotting is at the bottom axis. In the multivariate case, points are 
    plotted in the form of a scatterplot or matrix of scatterplots; this
    can be regulated by argument \code{data.par}.}
  
  \item{data.par}{an optional list of graphical parameters used for plotting 
    \code{data} in the multivariate case, when \code{data} is not \code{NULL}. 
    Recognized parameters are: \code{col}, \code{pch}, \code{cex}.
    If missing, the analogous components of \code{par()} are used. }
    
  \item{gap}{a numeric value which regulates the gap between panels of a
    multivariate plot when \code{d>2}.}
  
  \item{\dots}{additional graphical parameters}
  
}

\section{Details}{
  For univariate density plots, \code{probs} are used to compute quantiles
  from the appropriate distribution, and these are superimposed to the plot of
  the density function, unless \code{probs=NULL}. In the multivariate case,
  each bivariate plot  is constructed as a collection of contour curves,
  one curve for each probability level; consequently, \code{probs} cannot be 
  missing or \code{NULL}. The level of the density contour lines are chosen 
  so that each curve circumscribes a region with the  quoted probability, 
  to a good degree of approssimation; for additional information, see 
  Azzalini and Capitanio (2014), specifically Complement 5.2 and p.179, 
  and references therein.
}

\value{an invisible list. In the univariate case the list has three components:
the input object representing the distribution and two numeric vectors with 
the coordinates of the plotted density values. 
In the multivariate case, the first element of the list is the input object 
representing the distribution and all subsequent list elements are lists with 
components of the panels comprising the matrix plot; 
the elements of these sub-lists are: 
the vectors of \code{x} and \code{y} coordinates, the names of the variables, 
the density values at the \code{(x,y)} points, a vector of the density levels 
of the curves appearing in each panel plot, with the corresponding approximate
probability content as a vector attribute.}

\author{Adelchi Azzalini} 
 
\references{
 Azzalini, A. with the collaboration of Capitanio, A. (2014). 
 \emph{The Skew-Normal and Related Families}. 
 Cambridge University Press, IMS Monographs series.    
}

\seealso{\code{\link{makeSECdistr}}, \code{\link{summary.SECdistr}}, 
  \code{\link{dp2cp}}}

\section{Methods}{
\describe{

% \item{\code{signature(x = "ANY", y = "ANY")}}{Generic function: see
%      \code{\link[graphics]{plot}}.}


\item{\code{signature(x = "SECdistrUv")}}{Plot an object \code{x}
  of class \code{SECdistrUv}.}
 

\item{\code{signature(x = "SECdistrMv")}}{Plot an object \code{x}
  of class \code{SECdistrMv}.}
  
}}

\examples{
# d=1
f1 <- makeSECdistr(dp=c(3,2,5), family="SC", name="Univariate Skew-Cauchy")
plot(f1)
plot(f1, range=c(-3,40), probs=NULL, col=4)
#    
# d=2
Omega2 <- matrix(c(3, -3, -3, 5), 2, 2)    
f2 <-  makeSECdistr(dp=list(c(10,30), Omega=Omega2, alpha=c(-3, 5)),  
                family="sn", name="SN-2d", compNames=c("x1","x2"))
plot(f2)    
x2 <- rmsn(100, dp=slot(f2,"dp"))
plot(f2, main="Distribution 'f2'", probs=c(0.5,0.9), cex.main=1.5, col=2,
  cex=0.8, compLabs=c(expression(x[1]), expression(log(z[2]-beta^{1/3}))), 
  data=x2, data.par=list(col=4, cex=0.6, pch=5))

}
    
\keyword{methods}
\keyword{hplot}
