\name{factor.model.stat}
\alias{factor.model.stat}
\title{ Estimate Variance Matrix via Statistical Factors }
\description{
Creates a variance matrix based on the principal components of the
variables that have no missing values.
}
\usage{
factor.model.stat(x, weights = seq(0.5, 1.5, length.out = nobs), 
	output = "full", center = TRUE, frac.var = 0.5, iter.max = 1, 
	nfac.miss = 1, full.min = 20, reg.min = 40, sd.min = 20, 
	quan.sd = 0.9, tol = 0.001, zero.load = FALSE, 
	range.factors = c(0, Inf), constant.returns.okay = FALSE, 
	specific.floor = 0.1, floor.type = "quantile", verbose=2)
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
	a vector of observation weights, or \code{NULL}.

	Equal weights can be specified with \code{NULL} or with a 
	single positive number.

	Otherwise, the length must be equal to either 
	the original number of rows
	in \code{x} or the number of rows in \code{x} minus the number
	of rows that contain all missing values.
	}
  \item{output}{ 
	a character string indicating the form of the result.
	It must partially match one of: \code{"full"}, \code{"systematic"},
	\code{"specific"} or \code{"factor"}.
	}
  \item{center}{ 
	either a logical value or a numeric vector with length equal to
	the number of columns in \code{x}.
	If \code{center} is \code{TRUE}, then the mean of each column is 
	used as the center.
	If \code{center} is \code{FALSE}, then the center for each variable
	is taken to be zero.
	}
  \item{frac.var}{ 
	a control on the number of factors to use -- the number of factors
	is chosen so that the factors account for (just over) \code{frac.var}
	of the total variability.
	}
  \item{iter.max}{ 
	the maximum number of times to iterate the search for 
	principal factors of the variables with complete data.
	}
  \item{nfac.miss}{ 
	a vector of integers giving the number of factors to use in regressions
	for variables with missing values.
	The number of factors used is equal to the i-th element of 
	\code{nfac.miss} where i is the number of missing values for the
	variable.
	Thus the values in the vector should be non-increasing.
	The last value is used when the number of missing values is greater
	than the length of \code{nfac.miss}.
	}
  \item{full.min}{ 
	an integer giving the minimum number of variables that must have
	complete data.
	}
  \item{reg.min}{ 
	the minimum number of non-missing values for a variable in order
	for a regression to be performed on the variable.
	}
  \item{sd.min}{ 
	the minimum number of non-missing values for a variable in order
	for the standard deviation to be estimated from the data.
	}
  \item{quan.sd}{ 
	the quantile of the standard deviations to use for the standard
	deviation of variables that do not have enough data for the standard
	deviation to be estimated.
	}
  \item{tol}{ 
	a number giving the tolerance for the principal factor convergence
	(using the assets with full data).
	If the maximum change in uniquenesses (in the correlation scale) is
	less than \code{tol} from one iteration to the next, then convergence
	is assumed and the iterations end.
	}
  \item{zero.load}{ 
	a logical value.
	If \code{TRUE}, then loadings for variables with missing values are
	zero except for those estimated by regression.
	If \code{FALSE}, then loadings for variables with missing values are
	the average loading for the factor (when they are not estimated
	by regression).
	}
  \item{range.factors}{ 
	a numeric vector that gives the maximum and minimum number of
	factors that are allowed to be used.
	}
  \item{constant.returns.okay}{ 
	a logical vector: if \code{TRUE}, then a column with all of its
	non-missing values equal does not cause an error.

	if the true variance is thought to be non-zero, then
	a better alternative to setting this to \code{TRUE} is to set all
	the values in the column of \code{x} to be \code{NA}.
	}
  \item{specific.floor}{ 
	a number indicating how much uniquenesses should be adjusted upwards.
	The meaning of this number depends on the value of the
	\code{floor.type} argument.
	}
  \item{floor.type}{ 
	a character string that partially matches one of:
	\code{"quantile"} or \code{"fraction"}.
	
	If the value is \code{"quantile"}, then all uniquenesses are
	made to be at least as big as the \code{specific.floor} quantile
	of the uniquenesses.

	If the value is \code{"fraction"}, then all uniqueness are made
	to be at least \code{specific.floor}.
	}
  \item{verbose}{ 
	a number indicating the level of warning messages desired.
	This currently controls warnings:

	If at least 1, then a warning will be issued if all the values
	in \code{x} are non-negative.  
	In finance this is an indication that prices rather than returns
	are input (an easy mistake to make).

	If at least 1, then a warning will be issued if there are any
	assets with constant returns (unless \code{constant.returns.okay}
	is \code{FALSE} in which case an error is thrown).

	If at least 2, then a warning will be issued if there are any
	specific variances that are adjusted from being negative.
	}
}
\value{
if \code{output} is \code{"full"}, then a variance matrix with dimensions
equal to the number of columns in the input \code{x}.
This has two additional attributes: \code{number.of.factors} that says
how many factors are used in the model, and \code{timestamp} that gives
the date and time that the object was created.

if \code{output} is \code{"systematic"}, then a matrix with dimensions
equal to the number of columns in the input \code{x} that contains
the systematic portion of the variance matrix.

if \code{output} is \code{"specific"}, then a diagonal matrix with dimensions
equal to the number of columns in the input \code{x} that contains
the specific variance portion of the variance matrix.
The full variance matrix is the sum of the systematic and specific matrices.

If \code{output} is \code{"factor"}, then an object of class
\code{"statfacmodBurSt"} which is a list with components:
  \item{loadings }{
	a matrix of the loadings for the correlation matrix.
	}
  \item{uniquenesses }{
	the uniquenesses for the correlation matrix.
	That is, the proportion of the variance that is not explained
	by the factors.

	Note that if there are uniquenesses that have been modified
	via the \code{specific.floor} argument, then the actual proportion
	is the stated proportion divided by one plus the modification.
	}
  \item{sdev }{
	the standard deviations for the variables.

	Note that if there are uniquenesses that have been modified
	via the \code{specific.floor} argument, then the corresponding
	standard deviations in \code{sdev} are smaller than the actual
	standard deviations in the answer.
	}
  \item{constant.names }{
	A character vector giving the names of the variables that are
	constant (if any).
	}
  \item{cumulative.variance.fraction }{
	numeric vector giving the cumulative fraction of the variance 
	explained by (all) the factors.
	}
  \item{timestamp }{
	character string giving the date and time the calculation
	was completed.
	}
  \item{call }{
	an image of the call that created the object.
	}
}

\section{Details }{
Observations that are missing on all variables are deleted. 
Then a principal components factor model is estimated 
with the variables that have complete data.

For variables that have missing values, the standard deviation is estimated
when there are enough obeservations otherwise a given quantile of the 
standard deviations of the other assets is used as the estimate.
The loadings for these variables are set to be either 
the average loading for the variables with no missing data, or zero. 
The loadings for the most important factors are modified by performing a 
regression with the non-missing data for each variable (if there is enough
data to do the regression).

The treatment of variables with missing values can be quite important.
You may well benefit from specializing how missing values are handled to
your particular problem.
To do this, set the output to \code{"factor"} -- then you can modify the
loadings (and per force the uniquenesses), and the standard deviations to
fit your situation.
This may include taking sectors and countries into account, for example.

The default settings for missing value treatment are suitable for creating
a variance matrix for long-only portfolio optimization -- high volatility
and average correlation.
Take note that the proper treatment of missing values is HIGHLY dependent on
the use to which the variance matrix is to be put.

OBSERVATION WEIGHTS.
Time weights are quite helpful for estimating variances from returns.
The default weighting seems to perform reasonably well over a range of
situations.

FACTOR MODEL TO FULL MODEL.
This class of object has a method for \code{fitted} which
returns the variance matrix corresponding to the factor model representation.
}

\section{Warning }{
The default value for \code{weights} assumes that the last row is the most
recent observation and the first observation is the most ancient observation.
}

\section{Research Issues }{
The method of handling missing values used in the function has not
been well studied.
It seems not to be the worst approach, but undoubtedly can be improved.

The default method of boosting the result away from singularity
is completely unstudied.
For optimization it is wise to move away from singularity, just how to do
that best seems like a research question.
}

\section{Revision }{
This help was last revised 2014 March 09.
}

\author{ Burns Statistics }

\seealso{ \code{\link{fitted.statfacmodBurSt}}, 
	\code{\link{var.shrink.eqcor}}, \code{cov.wt},
	\code{\link{slideWeight}}.
}

\examples{
\dontrun{
varian1 <- factor.model.stat(retmat)

varfac <- factor.model.stat(retmat, nfac=0, zero=TRUE, output="fact")

varian2 <- fitted(varfac) # get matrix from factor model

varian3 <- factor.model.stat(retmat, nfac=rep(c(5,3,1), c(20,40,1)))
}
}
\keyword{multivariate}
