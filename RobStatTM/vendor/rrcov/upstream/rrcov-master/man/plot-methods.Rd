\name{plot-methods}
\docType{methods}
\concept{High breakdown point}
\alias{plot-methods}
\alias{plot,CovClassic-method}
\alias{plot,CovClassic,missing-method}
\alias{plot,CovRobust-method}
\alias{plot,CovRobust,missing-method}
\title{Methods for Function 'plot' in Package 'rrcov'}
\description{
 Shows the Mahalanobis distances based on robust and/or classical estimates 
of the location and the covariance matrix in different plots. 
The following plots are available:

 - index plot of the robust and mahalanobis distances
 
 - distance-distance plot
 
 - Chisquare QQ-plot of the robust and mahalanobis distances
 
 - plot of the tolerance ellipses (robust and classic)
 
 - Scree plot - Eigenvalues comparison plot
 
}

\usage{
\S4method{plot}{CovClassic}(x, which = c("all","distance","qqchi2","tolellipse","screeplot"), 
        ask=(which=="all" && dev.interactive()), 
        cutoff, id.n, tol=1e-7, \dots)
\S4method{plot}{CovRobust}(x, which = c("all","dd","distance","qqchi2","tolellipse","screeplot"), 
        classic=FALSE, ask=(which=="all" && dev.interactive()), 
        cutoff, id.n, tol=1e-7, \dots)
}

\arguments{
 \item{x}{an object of class \code{"Cov"}  or \code{"CovRobust"}}
  \item{which}{ Which plot to show? See Details for description of the options. Default is \code{which}="all". }.
  \item{classic}{ whether to plot the classical distances too. Default is \code{classic}=FALSE. }.
  \item{ask}{ logical; if 'TRUE', the user is \emph{ask}ed before each plot, see 'par(ask=.)'. 
Default is \code{ask = which=="all" && dev.interactive()}.  }
  \item{cutoff}{ The cutoff value for the distances.  }
  \item{id.n}{ Number of observations to identify by a label. If not supplied, the number of observations with distance larger than \code{cutoff} is used.  }
  \item{tol}{ tolerance to be used for computing the inverse see 'solve'. Default is \code{tol = 10e-7}}
  \item{... }{ other parameters to be passed through to plotting functions. }
}
\section{Methods}{\describe{

\item{x = "Cov", y = "missing"}{Plot mahalanobis distances for \code{x}.}
\item{x = "CovRobust", y = "missing"}{Plot robust and classical mahalanobis distances for \code{x}.}
}}
\examples{
data(hbk)
hbk.x <- data.matrix(hbk[, 1:3])
cv <- CovClassic(hbk.x)
plot(cv)
rcv <- CovMest(hbk.x)
plot(rcv)
}
\keyword{methods}
