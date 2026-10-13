\name{alpha.proxy}
\alias{alpha.proxy}
\title{ Compute and Plot Alpha Proxy }
\description{
Computes (and possibly generates a contour plot of) the alpha proxy -- a
measure of the effect that volatility and correlation have on the
utility of the investor.
As the name might suggest, it is in the same units as alpha (that is,
expected returns).
}
\usage{
alpha.proxy(weight = 0.2, vol.man = 0.2, vol.bench = 0.2, 
	vol.other = 0.2, cor.man = 0.2, cor.bench = 0.2, 
	plot.it = TRUE, transpose = FALSE, ...)
}
\arguments{
  \item{weight}{ 
	a number or vector of the fraction of the value of the investor's
	entire portfolio that is given to the manager.
	}
  \item{vol.man}{ 
	a number or vector giving the volatility of the manager's portfolio.
	}
  \item{vol.bench}{ 
	a number or vector giving the volatility of the benchmark.
	}
  \item{vol.other}{ 
	a number or vector giving the volatility of the rest of the portfolio.
	}
  \item{cor.man}{ 
	a number or vector giving the correlation between the manager's
	portfolio and the rest of the investor's portfolio.
	}
  \item{cor.bench}{ 
	a number or vector giving the correlation between the benchmark
	and the rest of the investor's portfolio.
	}
  \item{plot.it}{ 
	logical value.
	If \code{TRUE} and two of the arguments above have length greater
	than 1, then a contour plot is created.
	}
  \item{transpose}{ 
	logical value.
	If \code{TRUE}, then the variables on the axes of the contour
	plot are switched.
	}
  \item{\dots}{ 
	additional arguments to \code{filled.contour} may be given.
	}
}
\value{
a vector of the alpha proxies (in basis points)
if less than 2 of the first 6 arguments have length more than 1.

Otherwise it is a list (returned invisibly if \code{plot.it} is \code{TRUE})
with the following components:
  \item{x }{
	one of the vectors of inputs.
	}
  \item{y }{
	the other vector of inputs.
	}
  \item{z }{
	a matrix of the computed alpha proxies where rows correspond to
	the values in \code{x} and the columns correspond to \code{y}.
	The values are given in basis points.
	}
  \item{call }{
	a character string of the image of the command that created
	the object.
	}
}

\section{Side Effects }{
if \code{plot.it} is \code{TRUE}, then a contour plot is created.

An error occurs if more than 2 of the first 6 arguments have length
greater than 1.
}
\section{Details }{
The first 6 arguments are the variables that determine the alpha proxy.
The investor is faced with the decision of hiring the manager for some
portion (the weight) of the portfolio rather than using a replication
of the benchmark.
If the alpha proxy is positive, then the volatility and correlation of
the manager's portfolio is improving whatever outperformance the manager
may have (or is offsetting the losses).
}
\references{ 
Burns, Patrick (2003). "Portfolio Sharpening". Working Paper,
Burns Statistics http://www.burns-stat.com/.
	}
\section{Revision }{
This help was last revised 2010 January 05.
}

\seealso{ 
\code{\link{partial.rainbow}}.
}
\examples{
# return vector of alpha proxies
alpha.proxy(weigh=.05, vol.man=.17, cor.man=seq(0, .2, len=21))

# create a contour plot
alpha.proxy(weigh=.05, vol.man=seq(.15, .25, len=20), 
	cor.man=seq(0, .2, len=21))

# commands used to create figures in the paper
alpha.proxy(vol.man=seq(.15, .25, len=50), weight=seq(.01, .7, leng=60), 
        color.palette=partial.rainbow(start=0, end=.32))

alpha.proxy(cor.man=seq(0, .3, len=50), weight=seq(.01, .7, leng=60), 
        color.palette=partial.rainbow(start=.07))

alpha.proxy(cor.man=seq(0, .3, len=50), vol.man=seq(.15, .25, leng=60), 
        color.palette=partial.rainbow(start=0))
}
\keyword{dplot}
