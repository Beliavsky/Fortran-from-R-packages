%  file sn/man/st.prelimFit.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013-2020 Adelchi Azzalini
%
\name{st.prelimFit}
\alias{st.prelimFit}
\alias{mst.prelimFit}
\title{
Compute preliminary estimates for a linear model with ST-distributed error term
}
\description{
For a univariate or multivariate linear model where the error term
is assumed to have skew-\emph{t} (ST) distribution and the location parameter
is a linear function of a set of explanatory values, the functions compute
preliminary estimates to be used as initial values for a subsequent
maximization of the likelihood function.
These functions are mainly intended for internal package use.}

\usage{
st.prelimFit(x, y, w, quick = TRUE, verbose = 0, max.nu = 30, SN=FALSE)
mst.prelimFit(x, y, w, quick = TRUE, verbose = 0, max.nu = 30, SN=FALSE)
}

\arguments{
  \item{x}{design matrix of numeric values. It may be missing;
          if present, the first column must contain all 1's.}
  \item{y}{vector of observations of length \code{n},
     or a matrix with \code{n} rows.}
  \item{w}{a vector of non-negative integer weights of length \code{n}; 
    if missing, a vector of all 1's is generated.}
  \item{quick}{logical value which regulates whether a very quick estimate 
    is produced (default value \code{TRUE}); 
    see \sQuote{Details} for additional information.}
  \item{verbose}{an integer value which regulates the amount of messages
    printed out; default value is 0.}
  \item{max.nu}{threshold for the estimated degrees of freedom}
  \item{SN}{logical value (default value: \code{FALSE}); 
    if \code{TRUE}, a \acronym{SN} distribution is assumed.} 
}

\details{
The underlying methodology is the one presented by Azzalini and Salehi (2020).
In its essence, it is based on the selection of parameter values achieving the 
best matching between certain quantile-based summaries of the ST distribution
and the corresponding empirical quantities for the sample or, in the presence 
of explanatory variables, the same quantities computed from the residuals 
after fitting a median regression.

Argument \code{quick} selects whether the above-described matching is performed
in a quick or in an accurate way. Since the output values of this function are 
intended to be only initial values for subsequent likelihood maximization,
this explains the default option \code{quick=TRUE}. Other possible values
are \code{FALSE} and \code{NULL}; the latter simply sets \code{alpha=0}
and \code{nu=10}.

Since the methodology hinges on some selected sample quantiles, 
it can occasionally be spoiled by poor behaviour of these basic quantiles,
especially for small or moderate sample sizes.
The more visible effect of such situation is a very large value of the
estimated degrees of freedom, which then hampers rather than help a
subsequent likelihood maximization. It is therefore appropriate to set
an upper limit \code{max.nu} to this component.

Argument \code{x} may be missing. In this case, a one-column matrix with
all elements 1 is created.
}

\value{
A call to \code{st.prelimFit} returns a list with these components:
\item{dp}{a vector of estimates in the DP parameterization}
\item{residuals}{a vector of residual values}
\item{logLik}{the corresponding log-likelihood value}

A call to \code{mst.prelimFit} returns a list with these components:
\item{dp}{a list with the estimates in the DP parameterization}
\item{shrink.steps}{the number of shrinking steps applied to the original 
     estimate of the scale matrix to obtain an admissible matrix}
\item{dp.matrix}{a numeric matrix formed by the component-wise DP estimates}     
\item{logLik}{the corresponding log-likelihood value}
}
\references{
Azzalini, A. and Salehi, M. (2020).
Some computational aspects of maximum likelihood estimation of the 
skew-\emph{t} distribution.
In: \emph{Computational and Methodological Statistics and Biostatistics},
edited by  Andriëtte Bekker, Ding-Geng Chen and Johannes T. Ferreira.
Springer. DOI: 10.1007/978-3-030-42196-0
% https://www.springer.com/gp/book/9783030421953#aboutBook
}
\author{
Adelchi Azzalini}

\note{These functions are mainly intended to be called by \code{\link{selm}},
but they could be of interest for people developing their own procedures.}

\seealso{\code{\link[sn]{selm}} and  either \code{\link{dst}} or 
\code{\link{dmst}} for explanation of the DP parameters}

\examples{
data(barolo)
attach(barolo)  
A75 <- (reseller=="A" & volume==75)  
log.price <- log(price[A75], 10)
prelimFit <- st.prelimFit(y=log.price)
detach(barolo)
# 
data(ais)
attach(ais)
prelim32 <- mst.prelimFit(y=cbind(BMI, LBM), x=cbind(1, Ht, Wt))
detach(ais)
}
\keyword{robust} 
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
