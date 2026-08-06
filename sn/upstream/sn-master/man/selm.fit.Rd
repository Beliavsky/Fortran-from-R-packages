%  file sn/man/selm.fit.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013 Adelchi Azzalini
%---------------------
\name{selm.fit}
\alias{selm.fit}
\alias{sn.mple}
\alias{st.mple}
\alias{msn.mle}
\alias{msn.mple}
\alias{mst.mple}
 
\title{Fitting functions for \code{selm} models}

\description{A call to \code{selm} activates a call to \code{selm.fit} and
  from here to some other function which actually performs the parameter
  search, among those listed below. These lower-level functions can be 
  called directly for increased efficiency, at the expense of some more 
  programming effort and lack of methods for the returned object.}

\usage{ 
selm.fit(x, y, family = "SN", start = NULL, w, fixed.param = list(), 
   offset = NULL, selm.control=list())

sn.mple(x, y, cp = NULL, w, penalty = NULL, trace = FALSE, opt.method =
   c("nlminb",  "Nelder-Mead", "BFGS", "CG", "SANN"), control = list()) 

st.mple(x, y, dp = NULL, w, fixed.nu = NULL, symmetr = FALSE, penalty = NULL, 
   trace = FALSE, opt.method = c("nlminb", "Nelder-Mead", "BFGS", "CG", "SANN"), 
   control = list()) 

msn.mle(x, y, start = NULL, w, trace = FALSE, opt.method = c("nlminb", 
   "Nelder-Mead", "BFGS", "CG", "SANN"), control = list())

msn.mple(x, y, start = NULL, w, trace = FALSE, penalty = NULL, 
   opt.method = c("nlminb", "Nelder-Mead", "BFGS", "CG", "SANN"), 
   control = list()) 

mst.mple(x, y, start = NULL, w, fixed.nu = NULL, symmetr=FALSE, 
   penalty = NULL, trace = FALSE, 
   opt.method = c("nlminb", "Nelder-Mead", "BFGS", "CG", "SANN"), 
   control = list()) 
}

\arguments{
  \item{x}{a full-rank design matrix with the first column of all 1's.}

  \item{y}{a vector or a matrix of response values such that
     \code{NROW(y)=nrow(x)}.}

  \item{family}{a character string which selects the parametric family of
    distributions assumed for the error term of the regression model. 
    It must one of \code{"SN"} (default), \code{"ST"} or \code{"SC"}, which 
    correspond to the skew-normal, the skew-\emph{t} and the skew-Cauchy 
    family, respectively.
    See \code{\link{makeSECdistr}} for more information on these families and
    the skew-elliptically contoured (\acronym{SEC}) distributions; notice that
    family \code{"ESN"} is not allowed here.}

  \item{start, dp, cp}{a vector or a list of initial parameter values,
    depeding whether \code{y} is a vector or a matrix. It is assumed that 
    \code{cp} is given in the \acronym{CP} parameterization, \code{dp} and
    \code{start} in the \acronym{DP} parameterization. 
    For \code{st.mple} and \code{mst.mple}, see also the paragraph about  
    \code{start} in the documentation \sQuote{Details} of \code{selm}.
    }

  \item{w}{a vector of non-negative integer weights of length equal to 
    \code{NROW(y)}; if missing, a vector of all 1's is generated.}
    
  \item{fixed.param}{a list of assignments of parameter values to be kept
   fixed during the optimization process. Currently, there is only one such
   option, namely \code{fixed.param=list(nu='value')}, to fix the degrees
   of freedom at the named \code{'value'} when \code{family="ST"}, for instance
   \code{list(nu=3)}. Setting \code{fixed.param=list(nu=1)} is equivalent to
   select \code{family="SC"}.}
   
  \item{penalty}{an optional character string with the name of the penalty
   function of the log-likelihood; 
   default value \code{NULL} corresponds to no penalty.}

  \item{offset}{this can be used to specify an \emph{a priori} known
   component to be included in the linear predictor during fitting.  This
   should be \code{NULL} or a numeric vector of length equal to the number of
   cases.  One or more \code{\link{offset}} terms can be included in the
   formula instead or as well, and if more than one are specified their sum is
   used.}  %See \code{\link[stats]{model.offset}.}
 
  \item{trace}{a logical value which regulates printing of successive calls 
   to the target function; default value is \code{FALSE} which suppresses 
   printing.}
   
  \item{fixed.nu}{a positive value to keep fixed the parameter \code{nu} 
   of the \acronym{ST} distribution in the optimization process; with default
   value \code{NULL}, \code{nu} is estimated like the other parameters.}
   
  \item{symmetr}{a logical flag indicating whether a contraint of symmetry is
   imposed on the slant parameter; default is \code{symmetr=FALSE}.}  

  \item{opt.method}{a character string which selects the optimization method
    within the set \code{c("nlminb", "Nelder-Mead", "BFGS", "CG", "SANN")};
    the last four of these are \code{"methods"} of function \code{optim}.}

  \item{selm.control}{a list whose components regulate the working of 
    \code{selm.fit}; see \sQuote{Details} for their description;}

  \item{control}{a list of control items passed to the optimization function.}
}

\details{
A call to \code{selm} produces a call to \code{selm.fit} which
selects the appropriate function among \code{sn.mple}, \code{st.mple},
\code{msn.mle}, \code{msn.mple}, \code{mst.mple}, depending on the
arguments of the calling statement. In the adopted scheme for function names, 
\code{msn} refers to a  multivariate skew-normal distribution and  
\code{mst} refers to a multivariate skew-\eqn{t} distribution, while
\code{mle} and \code{mple} refers to maximum likelihood and maximum
penalized likelihood estimation, respectively.
Of these functions, \code{sn.mple} works in \acronym{CP} space; the others
in the \acronym{DP} space. In all cases, a correspondig mapping to the 
alternative parameter space is performed before exiting \code{selm.fit},
in addition to the selected parameter set.

The components of \code{selm.control} are as follows:
  \itemize{
    \item \code{method}: the estimation method, \code{"MLE"} or \code{"MPLE"}.
    \item \code{penalty}: a string with the name of the penalty function.
    \item \code{info.type}: a string with the name of the information matrix,
        \code{"observed"} or \code{"expected"}; currently fixed at "observed".         
    \item \code{opt.method}: a character string which selects the optimization
       method. 
    \item \code{opt.control}: a list of control parameters of \code{opt.method}.
    }
    
Function \code{msn.mle}, for \acronym{MLE} estimation of linear models with
\acronym{SN} errors, is unchanged from version 0.4-x of the package. 
Function \code{msn.mple} is similar to \code{msn.mle} but allows to introduce
a penalization of the log-likelihood; when \code{penalty=NULL}, a call to 
\code{msn.mle} is more efficient.
Functions \code{sn.mple} and \code{mst.mple} work like \code{sn.mle} and
\code{mst.mle} in version 0.4-x if the argument \code{penalty} is not 
set or it is set to \code{NULL}, except that \code{mst.mple} does not
handle a univariate response (use \code{st.mple} for that).
}

\value{A list whose specific components depend on the named function.
Typical components are:
 \item{call}{the calling statement}
 \item{dp}{vector or list of estimated \acronym{DP} parameters}
 \item{cp}{vector or list of estimated \acronym{CP} parameters}
 \item{logL}{the maximized (penalized) log-likelihood}
 \item{aux}{a list with auxiliary output values, depending on the function}
 \item{opt.method}{a list produced by the numerical \code{opt.method}}
}

\section{Background}{
Computational aspects of maximum likelihood estimation for univariate 
\acronym{SN} distributions are discussed in Section 3.1.7 of Azzalini and
Capitanio (2014). The working of \code{sn.mple} follows these lines; 
maximization is performed in the \acronym{CP} space. All other functions
operate on the \acronym{DP} space.
    
The technique underlying \code{msn.mle} is based on a partial analytical 
maximization, leading  implicitly to a form of profile log-likelihood.
This scheme is formulated in detail in Section 6.1 of Azzalini and Capitanio 
(1999) and summarized in Section 5.2.1 of Azzalini and Capitanio (2014). 
The same procedure is not feasible when one adopts \acronym{MPLE}; 
hence function \code{msn.mple} has to maximize over a larger parameter space.

When the SN family is fitted with the constraint \kbd{alpha=0}, this amounts
to adopt a classical linear model with Gaussian distributional assumption.
The corresponding \acronym{MLE}'s are the same as those produced by \code{lm}, 
except that the denominator the of the \acronym{MLE} variance (matrix) has the 
`uncorrected' form.
In the multivariate case, the covariance matrix of \acronym{MLE} is computed 
using expression (10) in Section 15.8 of Magnus and Neudecker (2007).

Maximization of the univariate \acronym{ST} log-likelihood is speeded-up 
by using the expressions of the gradient given by DiCiccio and Monti (2011),
reproduced with inessential variants in Section 4.3.3 of Azzalini and 
Capitanio (2014).

The working of \code{mst.mple} is based on a re-parameterization described 
in Section 5.1 of Azzalini and Capitanio (2003). The expressions of the 
corresponding log-likelihood derivatives are given in Appendix B of the full 
version of the paper.
}

\references{

  Azzalini, A. and Capitanio, A. (1999).
  Statistical applications of the multivariate skew normal distribution.
  \emph{J.Roy.Statist.Soc. B} \bold{61}, 579--602. 
  Full-length version available at \url{https://arXiv.org/abs/0911.2093}

  Azzalini, A. and Capitanio, A. (2003).
  Distributions generated by perturbation of symmetry with emphasis on
  a multivariate skew \emph{t} distribution.
  \emph{J.Roy. Statist. Soc. B} \bold{65}, 367--389.
  Full-length version available at \url{https://arXiv.org/abs/0911.2342}

  Azzalini, A. with the collaboration of Capitanio, A. (2014). 
  \emph{The Skew-Normal and Related Families}. 
  Cambridge University Press, IMS Monographs series.

  DiCiccio,  T. J. and Monti, A. C. (2011).  
  Inferential aspects of the skew \eqn{t}-distribution.
  \emph{Quaderni di Statistica} \bold{13}, 1--21.
  
  Magnus, J. R. and Neudecker, H. (2007).
  \emph{Matrix Differential Calculus with Applications in Statistics and 
  Econometrics}, third edition.  John Wiley & Sons.
}

\author{Adelchi Azzalini}

% \note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{selm}} for a comprehensive higher level fitting function,   
  \code{\link{Qpenalty}} for specification of a penalty function 
  }

\examples{ 
data(wines, package="sn")
X <- model.matrix(~ phenols + wine, data=wines)
fit <- msn.mle(x=X, y=cbind(wines$acidity, wines$alcohol), opt.method="BFGS")
fit <- st.mple(x=X, y = wines$acidity, fixed.nu=4, penalty="Qpenalty")
}
\keyword{regression}
\keyword{multivariate}
