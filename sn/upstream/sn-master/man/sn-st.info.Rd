%  file sn/man/sn-st.info.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013 Adelchi Azzalini
%---------------------
\name{sn-st.info}
\alias{sn.infoUv}
\alias{sn.infoMv}
\alias{st.infoUv}
\alias{st.infoMv}
\title{Expected and observed Fisher information for \acronym{SN} 
  and \acronym{ST} distributions}

\description{
  Computes Fisher information for parameters of simple sample having
  skew-normal (\acronym{SN}) or skew-\eqn{t} (\acronym{ST}) distribution or
  for a regression model with errors term having such distributions, in the
  \acronym{DP} and \acronym{CP} parametrizations.
}

\usage{
sn.infoUv(dp=NULL, cp=NULL, x=NULL, y, w, penalty=NULL,  norm2.tol=1e-06) 
    
sn.infoMv(dp, x=NULL, y, w, penalty=NULL, norm2.tol=1e-06, at.MLE=TRUE)

st.infoUv(dp = NULL, cp = NULL, x = NULL, y, w, fixed.nu = NULL, 
    symmetr = FALSE, penalty = NULL, norm2.tol = 1e-06) 

st.infoMv(dp, x = NULL, y, w, fixed.nu = NULL, symmetr = FALSE, 
    penalty = NULL, norm2.tol = 1e-06) 
}

\arguments{

\item{dp, cp}{direct or centred parameters, respectively; one of them
  can be a non-\code{NULL} argument. For the univariate \acronym{SN}
  distribution, \code{sn.infoUv} is to be used, and these arguments are
  vectors. In the multivariate case, \code{sn.infoMv} is to be used and these
  arguments are lists. See \code{\link{dp2cp}} for their description.}

\item{x}{an optional matrix which represents the design matrix of a 
  regression model}

\item{y}{a numeric vector (for \code{sn.infoUv} and \code{st.infoUv})
   or a matrix (for \code{sn.infoMv} and \code{st.infoMv}) representing the
   response.  In the \acronym{SN} case ( \code{sn.infoUv} and
   \code{sn.infoMv}), \code{y} can be missing, and in this case the expected
   information matrix is computed; otherwise the observed information is
   computed. In the \acronym{ST} case (\code{st.infoUv} and \code{st.infoMv}),
   \code{y} is a required argument, since only the observed information matrix
   for \acronym{ST} distributions is implemented. See \sQuote{Details} for
   additional information.}

\item{w}{an optional vector of weights (only meaningful for the observed 
   information, hence if \code{y} is missing); if missing, a vector of 1's is
   generated.}

\item{fixed.nu}{an optional numeric value which declares a fixed value of the
   degrees of freedom, \code{nu}. If not \code{NULL}, the information matrix
   has a dimension reduced by 1.}
   
\item{symmetr}{a logical flag which indicates whether a symmetry condition of 
   the distribution is being imposed; default is \code{symmetr=FALSE}.}   

\item{penalty}{a optional character string with the name of the penalty 
  function used in the call to \code{\link{selm}}; see this function for its
  description.}

\item{norm2.tol}{for the observed information case, the Mahalanobis squared
  distance of the score function from 0 is evaluated; if it exceeds
  \code{norm2.tol}, a warning message is issued, since the \sQuote{information
  matrix} so evaluated may be not positive-definite.  See \sQuote{Details} for
  additional information.}
  
\item{at.MLE}{a logical flag; if \code{at.MLE=TRUE} (default value), 
  it generates  a warning if the information matrix is not positive 
  definite or an error if the observed information matrix is not evaluated 
  at a maximum of the log-likelihood function.}
  
}

\value{
a list containing the following components:
\item{dp, cp}{one of the two arguments is the one supplied on input; 
  the other one matches the previous one in the alternative parametrization.}

\item{type}{the type of information matrix: "observed" or "expected".}

\item{info.dp, info.cp}{matrices of Fisher (observed or expected) 
  information in the two parametrizations.}
 
\item{asyvar.dp, asyvar.cp}{inverse matrices of Fisher information in the two
  parametrizations, when available; See \sQuote{Details} for additional
  information. }
 
\item{aux}{a list containing auxiliary elements, depending of the selected 
  function and the type of computation.}
}

\section{Details}{

In the univariate \acronym{SN} case, when \code{x} is not set, then a simple
random sample is assumed and a matrix \code{x} with a single column of all 
1's is constructed; in this case, the supplied vector \code{dp} or \code{cp}
must have length 3.  If \code{x} is set, then the supplied vector of parameters,
\code{dp} or \code{cp}, must have length \code{ncol(x)+2}.
In the multivariate case, a direct extension of this scheme applies.

If the observed information matrix is required, \code{dp} or \code{cp} should
represent the maximum likelihood estimates (MLE) for the given \code{y},
otherwise the information matrix may fail to be positive-definite and it
would be meaningless anyway. 
Therefore, the squared Mahalobis norm of the score vector is evaluated and compared with \code{norm2.tol}. 
If it exceeds this threshold, this is taken as an indication that the supplied
parameter list is not at the \acronym{MLE} and a warning message is issued.
The returned list still includes \code{info.dp} and \code{info.cp}, but in 
this case these represent merely the matrices of second derivatives;
\code{asyvar.dp} and \code{asyvar.cp} are set to \code{NULL}.

}

\section{Background}{
  The information matrix for the the univariate \acronym{SN} distribution in
  the two stated parameterizations in discussed in Sections 3.1.3--4 of 
  Azzalini and Capitanio (2014). For the multivariate distribution,   
  Section 5.2.2 of this monograph summarizes briefly the findings of 
  Arellano-Valle and Azzalini (2008).
 
  For \acronym{ST} distributions, only the observed information matrix is 
  provided, at the moment. Computation for the univariate case is based on 
  DiCiccio and Monti (2011). For the multivariate case, the score function is
  computed using an expression of Arellano-Valle (2010) followed by numerical
  differentiation.
}

\references{
  Arellano-Valle, R. B. (2010).
  The information matrix of the multivariate skew-\emph{t} distribution.
  \emph{Metron}, \bold{LXVIII}, 371--386.

  Arellano-Valle, R. B., and Azzalini, A. (2008).
  The centred parametrization for the multivariate skew-normal distribution.
  \emph{J. Multiv. Anal.} \bold{99}, 1362--1382.
  Corrigendum: \bold{100} (2009), 816.

  Azzalini, A. with the collaboration of Capitanio, A. (2014). 
  \emph{The Skew-Normal and Related Families}. 
  Cambridge University Press, IMS Monographs series.
  
  DiCiccio,  T. J. and Monti, A. C. (2011).  
  Inferential aspects of the skew \eqn{t}-distribution.
  \emph{Quaderni di Statistica} \bold{13}, 1--21.
}

\seealso{\code{\link{dsn}}, \code{\link{dmsn}}, \code{\link{dp2cp}}}

\examples{
infoE <- sn.infoUv(dp=c(0,1,5))           # expected information
set.seed(1); rnd <- rsn(100, dp=c(0, 1, 3))
fit <- selm(rnd~1, family="SN")
infoO <- sn.infoUv(cp=coef(fit), y=rnd)   # observed information
#
data(wines)
X <- model.matrix(~ pH + wine, data=wines)
fit <- sn.mple(x=X, y=wines$alcohol)
infoE <- sn.infoUv(cp=fit$cp, x=X)
infoO <- sn.infoUv(cp=fit$cp, x=X, y=wines$alcohol)
}
\keyword{distribution}
