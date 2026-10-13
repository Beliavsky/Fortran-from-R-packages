\name{var.shrink.eqcor}
\alias{var.shrink.eqcor}
\title{Ledoit-Wolf Shrinkage Variance Estimate }
\description{
Returns a variance matrix that shrinks towards the equal correlation 
matrix -- a Ledoit-Wolf estimate.
}
\usage{
var.shrink.eqcor(x, weights = seq(0.5, 1.5, length = nt), shrink = NULL, 
	center = TRUE, vol.shrink = 0, sd.min = 20, quan.sd = 0.9, 
	tol = 1e-4, compatible = FALSE, verbose=2)
}
\arguments{
  \item{x}{ 
	required.
	A numeric matrix.
	The rows are observations and the columns are the variables.
	In finance, this will be a matrix of returns where the rows are
	times and the columns are assets.
	For the default value of \code{weights} the most recent observation
	should be the last row.
	The number of columns may exceed the number of rows, and missing
	values are accepted.
	A column may even have all missing values.
}
  \item{weights}{ 
	a numeric vector giving the observation weights, or \code{NULL}.

	Equal weights can be specified with \code{NULL} or with a 
	single positive number.

	Otherwise, this must be a
	vector of non-negative numbers that is as long as the number
	of observations (rows of \code{x}).
}
  \item{shrink}{
	either \code{NULL} meaning the shrinkage will be estimated,
	or a single number meant to be in the range of zero to one.
}
  \item{center}{
	either a single logical value, or a numeric vector as long as
	the number of columns of \code{x}.
	This gives the expected value of each column. 
	A value of \code{TRUE} means the (weighted) mean of the columns
	of \code{x} is used.
	A value of \code{FALSE} means zero is used.
}
  \item{vol.shrink}{ 
	a number between zero and one (inclusive) that says how much to
	shrink the standard deviations towards the mean standard deviation.
}
  \item{sd.min}{ 
	a single integer giving the minimum number of observations needed
	in a column of \code{x} for the estimated standard deviation to be
	used.
}
  \item{quan.sd}{ 
	a single number in the range of zero to one (inclusive) that says
	what quantile of the (used) standard deviations to use as the
	estimate of standard deviation for columns without enough (or any)
	data.
}
  \item{tol}{ 
	a single number indicating how large the smallest eigenvalues
	of the answer should be.
	All eigenvalues are made to be no smaller than \code{tol} times the
	maximum eigenvalue.
}
  \item{compatible}{ 
	a logical value.
	If \code{TRUE}, then the variance is scaled by the number of
	observations rather than one less than the number of observations,
	and there will be a warning if \code{weights} are not all the same.
	This makes the results compatible with Ledoit-Wolf computations
	and with the \code{cov.shrink} function in package \code{tawny}.
	This also makes the variances smaller, which is not necessarily
	a good thing.
}
  \item{verbose}{ 
	a number indicating the level of warning messages desired.
	This currently controls only one warning:

	If at least 1, then a warning will be issued if all the values
	in \code{x} are non-negative.  
	In finance this is an indication that prices rather than returns
	are input (an easy mistake to make).

}
}
\value{
	an estimate of the variance matrix of \code{x}.
	The sample variance is shrunk towards equal correlation.

	This has two additional attributes:
  \item{shrink }{
	the estimated or input amount of shrinkage towards the equal
	correlation matrix from the sample variance.
	}
  \item{timedate }{
	the date and time at which the computation was done.
	}
}

\section{Details }{
Time weights are quite helpful for estimating variances from returns.
The default weighting seems to perform reasonably well over a range of
situations.
However, time weighting was not studied for this estimator.
}

\section{Warning }{
The default value for \code{weights} assumes that the last row is the most
recent observation and the first observation is the most ancient observation.
}

\section{Research Issues }{
The method of handling missing values used in the function has not
been studied (at all).

The method of boosting the result away from singularity
is completely unstudied.
For optimization it is wise to move away from singularity, just how to do
that best seems like a research question.
The method used boosts the smallest eigenvalues, it might be better
to increase the diagonal.
}
\references{ 
Olivier Ledoit and Michael Wolf (2004) "Honey, I shrunk the sample 
covariance matrix". The Journal of Portfolio Management, volume 30,
number 4.
}

\section{Revision }{
This help was last revised 2014 March 09.
}

\author{ Burns Statistics }

\seealso{ \code{\link{factor.model.stat}}, \code{cov.wt},
	\code{\link{slideWeight}}.
}
\examples{
\dontrun{
var1 <- var.shrink.eqcor(return.matrix)

var.unweighted <- var.shrink.eqcor(return.matrix, weights=1)
}
}
\keyword{multivariate}
