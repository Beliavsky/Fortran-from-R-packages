\name{CEoptim}
\alias{CEoptim}
\title{Cross-Entropy optimizer}
\usage{
	CEoptim(f, f.arg=NULL, maximize=FALSE, continuous=NULL, discrete=NULL,
	    N=100L, rho=0.1, iterThr=1e4L, noImproveThr=5, verbose=FALSE)}


\arguments{
	\item{f}{Function to be optimized. Can have continuous
  and discrete arguments}
\item{f.arg}{List of additional fixed arguments passed to function f.}
\item{maximize}{Logical value determining whether to maximize or minimize the objective function}

\item{continuous}{List of arguments for the continuous optimization part consisting of:}
\itemize{
\item \code{mean} Vector of initial means.
\item \code{sd} Vector of initial standard deviations.
\item \code{smoothMean} Smoothing parameter for the vector of means.  Default value 1 (no smoothing).
\item \code{smoothSd}  Smoothing parameter for the standard
deviations. Default value 1 (no smoothing).
\item \code{sdThr} Positive numeric convergence threshold.  Check whether  
 the maximum standard deviation is smaller than \code{sdThr}. Default value 0.001.
 \item \code{conMat} Coefficient matrix of linear constraint \code{conMat} \eqn{x \le} \code{conVec}.   
\item \code{conVec}  Value vector of linear constraint \code{conMat} \eqn{x \le} \code{conVec}.  
}

\item{discrete}{List of arguments for the discrete
optimization part, consisting of:}
\itemize{
\item \code{categories} Integer vector which defines the allowed
values of the categorical variables. The \code{i}th categorical
variable takes values in the set \{0,1,\ldots,\code{categories(i)}-1\}.
\item \code{probs} List of  initial probabilities for the
categorical variables. Defaults to equal (uniform) probabilities.
\item \code{smoothProb} Smoothing parameter for the probabilities of
the categorical sampling distribution. Default value 1 (no smoothing).
\item \code{ProbThr}  Positive numeric convergence
threshold. Check whether all probabilities 
in the categorical sampling distributions
deviate less than \code{ProbThr} from either 0 or 1. Default value 0.001.}

\item{N}{Integer representing the CE sample size.}

\item{rho}{Value between 0 and 1 representing the elite
proportion.}


\item{iterThr}{Termination threshold on the 
largest number of iterations.}

\item{noImproveThr}{Termination threshold on the largest
number of iterations during which no improvement of the best function
value is found.}

\item{verbose}{Logical value set for CE progress
output.}}

\value{\code{CEoptim} returns an object of class "CEoptim" which is a list with the following components.
	\itemize{
\item{\bold{optimum} }  Optimal value of \code{f}. 
\item{\bold{optimizer} } List of the location of the optimal value, consisting
of:
\itemize{
 \item{\bold{continuous}  } Continuous part of the optimizer. 
 \item{\bold{discrete}  }  Discrete part of the optimizer.
 }
\item{\bold{termination}  } List of termination information consisting of:
\itemize{
  \item{\bold{niter}  } Total number of iterations upon termination.
  \item{\bold{convergence}  } One of the following statements:
  \itemize{
	\item{\bold{\code{Not converged}}}, if the number of iterations reaches
	\code{iterThr};
	\item{\bold{\code{The optimum did not change for noImproveThr
	iterations}}}, if the best value has not improved for
\code{noImproveThr} iterations;
\item{\bold{\code{Variances converged}}}, otherwise.
	}
  }
  \item{\bold{states}   } List of intermediate results computed at each
iteration. 
 It consists of the iteration number (\code{iter}), the best overall
 value
 (\bold{\code{optimum}}) and the worst value of the elite
 samples, (\code{gammat}). The means (\code{mean}) and maximum standard deviations
 (\code{maxSd}) of the elite set are also included for continuous
 cases, and the maximum deviations (\code{maxProbs}) of the sampling
	probabilities to either 0 or 1 are included for discrete cases.

  
\item{\bold{states.probs}  } List of categorical sampling probabilities computed
at each iteration. Will only be returned for discrete and mixed cases.
}}
	
\description{
\code{CEopt} is an optimization function based on the Cross-Entropy method
}
\note{
 Although partial parameter passing is allowed outside lists, it is
  recommended that parameters names are specified in full. Parameters
  inside lists have to specified completely.

  Because \code{CEoptim} is a random function it is useful to (1)
  set the seed for the random number generator (for testing purposes),
  and (2)
  investigate the quality of the results by repeating
  the optimization a number of times.
  }

\author{Tim Benham, Qibin Duan, Dirk P. Kroese, Benoit Liquet}

\references{ Benham T.,  Duan Q., Kroese D.P., Liquet B. (2017) \bold{CEoptim}:
Cross-Entropy \R package for optimization. \emph{Journal of Statistical
  Software}, 76(8), 1-29.

Rubinstein R.Y. and Kroese D.P. (2004). \emph{The Cross-Entropy
Method}. Springer, New York.
}

\examples{
## Maximizing the Peaks Function


fun <- function(x){
return(3*(1-x[1])^2*exp(-x[1]^2 - (x[2]+1)^2)
	-10*(x[1]/5-x[1]^3 - x[2]^5)*exp(-x[1]^2 - x[2]^2)
	-1/3*exp(-(x[1]+1)^2 - x[2]^2))}

set.seed(1234)

mu0 <- c(-3,-3); sigma0 <- c(10,10)
 
res <- CEoptim(fun,continuous=list(mean=mu0, sd=sigma0), maximize=TRUE)

## To extract the Optimal value of fun
res$optimum
## To extract the location of the optimal value
res$optimizer$continuous
## print function gives the following default values
print(res)
}

