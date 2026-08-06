%  file sn/man/selm.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013-2017 Adelchi Azzalini
%---------------------
\name{selm}
\encoding{UTF-8}
\alias{selm}
\concept{regression}
\concept{skew-elliptical distribution}

\title{Fitting linear models with skew-elliptical error term}

\description{Function \code{selm} fits a \code{l}inear \code{m}odel
  with \code{s}kew-\code{e}lliptical error term. 
  The term \sQuote{skew-elliptical distribution} is an abbreviated equivalent 
  of skew-elliptically contoured (\acronym{SEC}) distribution.
  The function works for univariate and multivariate response variables.}

\usage{
selm(formula, family = "SN", data, weights, subset, na.action, 
  start = NULL, fixed.param = list(), method = "MLE",  penalty=NULL, 
  model = TRUE, x = FALSE, y = FALSE, contrasts = NULL, offset, ...)
}
 
\arguments{

  \item{formula}{an object of class \code{"\link[stats]{formula}"} 
   (or one that can be coerced to that class): a symbolic description of the
   model to be fitted, using the same syntax used for the similar parameter of
   e.g. \code{"\link[stats]{lm}"}, with the restriction that the constant
   term must not be removed from the linear predictor. 
   % The details of model specification are given under \sQuote{Details}.  
   }

  \item{family}{a character string which selects the parametric family
   of \acronym{SEC} type  assumed for the error term. It must be one of
   \code{"SN"} (default), \code{"ST"} or \code{"SC"}, which correspond to the
   skew-normal, the skew-\emph{t} and the skew-Cauchy family, respectively.
   See \code{\link{makeSECdistr}} for more information on these families and
   the set of \acronym{SEC} distributions; notice that the family \code{"ESN"} 
   listed there is not allowed here.}

  \item{data}{an optional data frame  containing the variables in
   the model.  If not found in \code{data}, the variables are taken from
   \code{environment(formula)}, typically the environment from which
   \code{selm} is called.}

  \item{weights}{a numeric vector of weights associated to  individual
   observations. Weights are supposed to represent frequencies, hence must be
   non-negative integers (not all 0) and \code{length(weights)} must equal the
   number of observations. If not assigned, a vector of all 1's is generated.}

  \item{subset}{an optional vector specifying a subset of observations
   to be used in the fitting process. It works like the same parameter
   in \code{\link[stats]{lm}}.}

  \item{na.action}{a function which indicates what should happen
   when the data contain \code{NA}s.  The default is set by the
   \code{na.action} setting of \code{\link[base]{options}}.  
   The \sQuote{factory-fresh} default is \code{\link{na.omit}}.  
   Another possible value is \code{NULL}, no action.  
   % Value \code{\link[stats]{na.exclude}} can be useful.
   }

  \item{start}{a vector (in the univariate case) or a list (in the  
   multivariate case) of initial \acronym{DP} values for searching the 
   parameter estimates. See \sQuote{Details} about a choice of
   \kbd{start} to be avoided.  If \code{start=NULL} (default), 
   initial values are  selected by the procedure.
   If \code{family="ST"}, an additional option exists; see \sQuote{Details}.}

  \item{fixed.param}{a list of assignments of parameter values which must
   be kept fixed in the estimation process. 
   Currently, there only two types of admissible constraint: one is to
   set \code{alpha=0} to impose a symmetry condition of the distribution; 
   the other is to set \code{nu=<value>}, to fix the degrees of freedom  
   at the named \code{<value>} when \code{family="ST"}, for instance
   \code{list(nu=3)}.  See \sQuote{Details} for additional information.
   }

  \item{method}{a character string which selects the estimation method to be
   used for fitting. Currently, two options exist: \code{"MLE"} (default) and
   \code{"MPLE"}, corresponding to standard maximum likelihood and maximum
   penalized likelihood estimation, respectively. See \sQuote{Details} for
   additional information.  }

  \item{penalty}{a character string which denotes the penalty function to be
   subtracted to the log-likelihood function, when \code{method="MPLE"}; if
   \code{penalty=NULL} (default), a pre-defined function is adopted. See
   \sQuote{Details} for a description of the default penalty function and for
   the expected format of alternative specifications.  When
   \code{method="MLE"}, no penalization is applied and this argument has no
   effect.}


  \item{model, x, y}{logicals.  If \code{TRUE}, the corresponding components
   of the fit are returned.}
 
  \item{contrasts}{an optional list. See the \code{contrasts.arg} of
   \code{\link[stats]{model.matrix.default}}.}
 
  \item{offset}{this can be used to specify an \emph{a priori} known
   component to be included in the linear predictor during fitting.  This
   should be \code{NULL} or a numeric vector of length equal to the number of
   cases.  One or more \code{\link{offset}} terms can be included in the
   formula instead or as well, and if more than one are specified their sum 
   is used. }
   
  \item{\dots}{optional control parameters, as follows.
   \itemize{

    \item \code{trace}: a logical value which indicates whether intermediate
        evaluations of the optimization process are printed (default:
        \code{FALSE}).
    \item \code{info.type}: a character string which indicates the type of
        Fisher information matrix; possible values are \code{"observed"}
        (default) and \code{"expected"}. Currently, \code{"expected"} is
        implemented only for the \acronym{SN} family.

    \item \code{opt.method}: a character string which selects the numerical
        optimization method, among the possible values 
        \code{"nlminb", "Nelder-Mead", "BFGS", "CG", "SANN"}. 
        If \code{opt.method="nlminb"} (default),
        function \code{\link[stats]{nlminb}} is called, 
        otherwise function \code{\link[stats]{optim}} is called with 
        \code{method} equal to \code{opt.method}.

    \item \code{opt.control}: a list of control parameters which is passed on
        either to \code{nlminb} or to \code{optim}, depending on the chosen
        \code{opt.method}.
    }
  }
}

\details{By default, \code{selm} fits the selected model by maximum
  likelihood estimation (\acronym{MLE}), making use of some numerical
  optimization method.  Maximization is performed in one
  parameterization, usually \acronym{DP}, and then the estimates are mapped to
  other parameter sets, \acronym{CP} and pseudo-\acronym{CP}; 
  see \code{\link{dp2cp}} for more information on parameterizations. 
  These parameter transformations are carried out trasparently to the user. 
  The observed information matrix is used to obtain the estimated variance 
  matrix of the \acronym{MLE}'s and from this the standard errors.  
  Background information on \acronym{MLE} in the context of \acronym{SEC} 
  distributions is provided by Azzalini and Capitanio (2014); 
  see specifically Chapter 3, Sections 4.3, 5.2,  6.2.5--6. For additional
  information, see the original research work referenced therein as well as
  the sources quoted below.
 
  Although the density functionof SEC distributions are expressed using
  \acronym{DP} parameter sets, the methods associated to the objects created
  by this function communicate, by default, their outcomes in the \acronym{CP}
  parameter set, or its variant form pseudo-\acronym{CP} when \acronym{CP}
  does not exist; the \sQuote{Note} at \code{\link{summary.selm}} explains why. 
  A more detailed discussion  is provided by Azzalini and Capitanio 
  (1999,  Section 5.2) and Arellano-Valle and  Azzalini (2008, Section 4), 
  for the univariate and the multivariate SN case, respectively; 
  an abriged account is available in Sections 3.1.4--6 and 5.2.3 of 
  Azzalini and Capitanio (2014). For the ST case, see Arellano-Valle 
  and  Azzalini (2013).
  
  There is a known open issue which affects computation of the information
  matrix of the multivariate skew-normal distribution when the slant
  parameter \eqn{\alpha} approaches the null vector; see p.149 of
  Azzalini and Capitanio (2014). Consequently, if a model with
  multivariate response is fitted with \code{family="SN"} and the estimate
  \code{alpha} of \eqn{\alpha} is at the origin or neary so, the
  information matrix and the standard errors are not computed and a
  warning message is issued. In this unusual circumstance, a simple
  work-around is to re-fit the model with \code{family="ST"}, which will
  work except in remote cases when (i) the estimated degrees of freedom
  \code{nu} diverge and (ii) still \code{alpha} remains at the origin.

  The optional argument \code{fixed.param=list(alpha=0)} imposes the
  constraint \eqn{\alpha=0} in the estimation process; in the multivariate 
  case, the expression is interpreted in the sense that all the components  
  of vector \eqn{\alpha} are zero, which implies symmetry of the
  error distribution, irrespectively of the parameterization 
  subsequently adopted for summaries and diagnostics.
  When this restriction is selected, the estimation method cannot be
  set to \code{"MPLE"}. Under the constraint \eqn{\alpha=0},
  if \code{family="SN"}, the model is  fitted similarly to \code{lm}, except
  that here \acronym{MLE} is used for estimation of the covariance matrix. 
  If \code{family="ST"} or \code{family="SC"}, a symmetric Student's \eqn{t} 
  or Cauchy distribution is adopted. 
  
  Under the constraint \eqn{\alpha=0}, the location parameter \eqn{\xi}
  coincides with the mode and the mean of the distribution, when the latter
  exists. In addition, when the covariance matrix of a \acronym{ST}
  distribution exists, it differs from \eqn{\Omega} only by a multiplicative
  factor. Consequently, the summaries of a model of this sort automatically
  adopt the \acronym{DP} parametrization.
  
  The other possible form of constraint allows to fix the degrees of
  freedom when \code{family="ST"}. The two constraints can be combined 
  writing, for instance,  \code{fixed.param=list(alpha=0, nu=6)}.
  The constraint \code{nu=1} is equivalent to select \code{family="SC"}.
  In practice, an expression of type \code{fixed.param=list(..)} can be
  abbreviated to \code{fixed=list(..)}.
  
  Argument \kbd{start} allows to set the initial values, with respect to the
  \acronym{DP} parameterization, of the numerical optimization. 
  However, there is a specific choice of start to be avoided.
  When \kbd{family="SN"}, do not set the shape parameter \kbd{alpha} 
  exactly at 0, as this would blow-up computation of the log-likelihood 
  gradient and the Hessian matrix. This is not due to a software bug, 
  but to a known peculiar behaviour of the log-likelihood  function at 
  that specific point. Therefore, in the univariate case for instance, 
  do not set e.g. \kbd{start=c(12, 21, 0)}, but set instead something
  like \kbd{start=c(12, 21, 0.01)}. 
  % Also, setting such an initial $\alpha=0$ or close to 0 is a questionable 
  % choice anyway: if one fits a model of this class, then  some asymmetry is
  % expected to be present and it is odd to start the search from a symmetry
  % condition.
  Recall that, if one needs to fit a model forcing 0 asymmetry, typically to
  compare two log-likelihood functions with/without asymmetry, then the option
  to use is \kbd{fixed.param=list(alpha=0)}.
  
  Since version 1.6.0, a new initialization procedure has been introduced
  for the case \kbd{family="ST"}, which adopts the method proposed by 
  Azzalini & Salehi (2020), implemented in functions \kbd{st.prelimFit} 
  and \kbd{mst.prelimFit}.
  Correspondingly, the \kbd{start} argument can now be of different type,
  namely a character with possible values \kbd{"M0"}, \kbd{"M2"} (detault in 
  the univariate case) and \kbd{"M3"} (detault in the multivariate case).
  The choice \kbd{"M0"} selects the older method, in use prior to version
  1.6.0. For more information, see Azzalini &  Salehi (2020).

  In some cases, especially for small sample size, the \acronym{MLE} occurs on
  the frontier of the parameter space, leading to \acronym{DP} estimates with
  \code{abs(alpha)=Inf} or to a similar situation in the multivariate case 
  or in an alternative parameterization. Such outcome is regared by many as
  unsatisfactory; surely it prevents using the observed information matrix to
  compute standard errors. This problem motivates the use of maximum penalized
  likelihood estimation (\acronym{MPLE}), where the regular log-likelihood
  function \eqn{\log~L}{log(L)} is penalized by subtracting an amount
  \eqn{Q}, say, increasingly large as \eqn{|\alpha|} increases. 
  Hence the function which is maximized at the optimization stage is now
  \eqn{\log\,L~-~Q}{log(L) - Q}.  If \code{method="MPLE"} and
  \code{penalty=NULL}, the default function \code{Qpenalty} is used,
  which implements the penalization:
     \deqn{Q(\alpha) = c_1 \log(1 + c_2 \alpha_*^2)}{%
           Q(\alpha)= c₁ log(1 + c₂ [\alpha*]²)}
  where \eqn{c_1}{c₁} and \eqn{c_2}{c₂} are positive constants, which
  depend on the degrees of freedom \code{nu} in the \code{ST} case,
      \deqn{\alpha_*^2 = \alpha^\top \bar\Omega \alpha}{%?
            [\alpha*]² = \alpha' cor(\Omega) \alpha}
  and \eqn{\bar\Omega}{cor(\Omega)} denotes the correlation matrix 
  associated to the scale matrix \code{Omega} described in connection with
  \code{\link{makeSECdistr}}. In the univariate case 
  \eqn{\bar\Omega=1}{cor(\Omega)=1},
  so that \eqn{\alpha_*^2=\alpha^2}{[\alpha*]²=\alpha²}. Further information 
  on \acronym{MPLE} and this choice of the penalty function is given in 
  Section 3.1.8 and p.111 of Azzalini and Capitanio (2014); for a more 
  detailed account, see Azzalini and Arellano-Valle (2013) and references  
  therein.

  It is possible to change the penalty function, to be declared via the 
  argument \code{penalty}. For instance, if the calling statement includes 
  \code{penalty="anotherQ"}, the user must have defined  

     \verb{    }\code{anotherQ <- function(alpha_etc, nu = NULL, der = 0)}

  with the following arguments.
  \itemize{
  \item \code{alpha_etc}: in the univariate case, a single value \code{alpha};
     in the multivariate case, a two-component list whose first component is
     the vector \code{alpha}, the second one is matrix equal to
     \code{cov2cor(Omega)}.
     % \eqn{\bar\Omega}{corOmega}.
  \item \code{nu}: degrees of freedom, only relevant if \code{family="ST"}.
  \item \code{der}: a numeric value which indicates the required order of
     derivation; if \code{der=0} (default value), only the penalty \code{Q}
      needs to be retuned by the function; 
      if \code{der=1}, \code{attr(Q, "der1")} must represent the
     first order derivative of \code{Q} with respect to \code{alpha}; if
     \code{der=2}, also \code{attr(Q, "der2")} must be assigned, containing
     the second derivative (only required in the univariate case).
    }
  This function must return a single numeric value, possibly with required
  attributes when is called with \code{der>1}.
  Since \pkg{sn} imports functions \code{\link[numDeriv]{grad}} and 
  \code{\link[numDeriv]{hessian}} from package \pkg{numDeriv}, one can rely 
  on them for numerical evaluation of the derivatives, if they are not 
  available in an explicit form.

  This penalization scheme allows to introduce a prior distribution 
  \eqn{\pi} for \eqn{\alpha} by setting \eqn{Q=-\log\pi}{Q=-log(\pi)}, 
  leading to a maximum \emph{a posteriori} estimate in the stated sense. 
  See \code{\link{Qpenalty}} for more information and an illustration.
  
  The actual computations are not performed within \code{selm} which only 
  sets-up ingredients for work of \code{\link{selm.fit}} and other functions
  further below this one.  See \code{\link{selm.fit}} for more information.
}

\value{an S4 object of class \code{selm} or \code{mselm}, depending on whether
  the response variable of the fitted model is univariate or multivariate;
  these objects are described in the \code{\linkS4class{selm} class}.
}

\references{
Arellano-Valle, R. B., and Azzalini, A. (2008).
 The centred parametrization for the multivariate skew-normal distribution.
 \emph{J. Multiv. Anal.} \bold{99}, 1362--1382.
  Corrigendum: \bold{100} (2009), 816.

Arellano-Valle, R. B., and Azzalini, A. (2013, available online 12 June 2011).
  The centred parametrization and related quantities for the  skew-\emph{t} 
  distribution.
  \emph{J. Multiv. Anal.} \bold{113}, 73--90. 

Azzalini, A. and Capitanio, A. (1999).
  Statistical applications of the multivariate skew normal distribution.
  \emph{J.Roy.Statist.Soc. B} \bold{61}, 579--602. 
  Full-length version available at \url{https://arXiv.org/abs/0911.2093}

Azzalini, A. and Arellano-Valle, R. B. (2013, available online 30 June 2012). 
  Maximum penalized likelihood estimation for skew-normal and skew-\emph{t} 
  distributions. \emph{J. Stat. Planning & Inference} \bold{143}, 419--433. 

Azzalini, A. with the collaboration of Capitanio, A. (2014). 
 \emph{The Skew-Normal and Related Families}. 
 Cambridge University Press, IMS Monographs series.

Azzalini, A. and Salehi, M. (2020).
  Some computational aspects of maximum likelihood estimation 
  of the skew-\emph{t} distribution. In \emph{Computational and Methodological 
  Statistics and Biostatistics}, edited by A. Bekker, Ding-Geng Chen and 
  Johannes T. Ferreira, pp.3-28.  Springer Nature Switzerland.

% Magnus and Neudecker
  
}

\author{Adelchi Azzalini}

\section{Cautionary notes}{
The first of these notes applies to the stage \emph{preceding} the
use of \kbd{selm} and related fitting procedures. Before fitting a model of
this sort, consider whether you have enough data for this task. 
In this respect, the passage below taken from p.63 of Azzalini 
and Capitanio (2014) is relevant.

\dQuote{Before entering technical aspects, it is advisable to underline 
a qualitative effect of working with a parametric family which effectively 
is regulated by moments up to the third order. 
The implication is that the traditional rule of thumb by which a sample 
size is small up to ‘about \eqn{n = 30}’, and then starts to become ‘large’,
while sensible for a normal population or other two-parameter distribution, 
is not really appropriate here. 
To give an indication of a new threshold is especially difficult, 
because the value of  \eqn{\alpha} also has a role here. 
Under this \emph{caveat}, numerical experience suggests that ‘about 
\eqn{n = 50}’ may be a more appropriate guideline in this context.}

The above passage referred to the univariate SN context. 
In the multivariate case, increase the sample size appropriately, 
especially so with the \acronym{ST} family.
This is not to say that one cannot attempt fitting these models 
with small or moderate sample size. However, one must be aware of the 
implications and not be surprised if problems appear.

The second cautionary note refers instead to the outcome of a call to 
\kbd{selm} and related function, or the lack of it.
The estimates are obtained by numerical optimization methods and, as
usual in similar cases, there is no guarantee that the maximum of the
objective function is achieved. Consideration of model simplicity
and of numerical experience indicate that models with \acronym{SN} error
terms generally produce more reliable results compared to those with 
the \acronym{ST} family. Take into account that models involving a 
traditional Student's \eqn{t} distribution with unknown degrees of freedom 
can already be problematic; the presence of the (multivariate) slant parameter
\eqn{\alpha} in the \acronym{ST} family cannot make things any simpler. 
Consequently, care must be exercised, especially so if one works with 
the (multivariate) \acronym{ST} family. 
Consider re-fitting a model with different starting values and, 
in the \acronym{ST} case, building the profile log-likelihood for a range 
of \eqn{\nu} values; function \code{\link{profile.selm}} can be useful here.

Details on the numerical optimization which has produced object \code{obj} 
can be extracted with \code{slot(obj, "opt.method")}; inspection of this
component can be useful in problematic cases.
# Be aware that  occasionally \code{optim} and \code{nlminb} declare successful
# completion of a regular minimization problem at a point where the Hessian 
# matrix is not positive-definite. 
}


\seealso{\itemize{

\item
  \code{\linkS4class{selm}-class} for classes \code{"selm"} and \code{"mselm"},
  \code{\link{summary.selm}} for summaries, \code{\link{plot.selm}} for plots,
   \code{\link{residuals.selm}} for residuals and fitted values

\item
  the generic functions \code{\link[=coef.selm]{coef}}, 
  \code{\link{logLik}}, \code{\link{vcov}}, 
  \code{\link[=profile.selm]{profile}}, 
  \code{\link[=confint.selm]{confint}}, 
  \code{\link[=predict.selm]{predict}}

\item
  the underlying function \code{\link{selm.fit}} and those further down

\item
  the selection of a penalty function of the log-likelihood, 
  such as \code{\link{Qpenalty}}
  
\item
  the function \code{\link{extractSECdistr}} to extract the \acronym{SEC}
  error distribution from an object returned by \code{selm} 
  
\item the broad underlying logic and a number of ingredients are like in 
   function \code{\link[stats]{lm}}  
}}

\examples{
data(ais)
m1 <- selm(log(Fe) ~ BMI + LBM, family="SN", data=ais)
print(m1)
summary(m1)
s <- summary(m1, "DP", cov=TRUE, cor=TRUE)
plot(m1)
plot(m1, param.type="DP")
logLik(m1)
coef(m1)
coef(m1, "DP")
var <- vcov(m1)
#
m1a <- selm(log(Fe) ~ BMI + LBM, family="SN", method="MPLE", data=ais)
m1b <- selm(log(Fe) ~ BMI + LBM, family="ST", fixed.param=list(nu=8), data=ais)
#
data(barolo)
attach(barolo)
A75 <- (reseller=="A" & volume==75)
logPrice <- log(price[A75],10) 
m <- selm(logPrice ~ 1, family="ST", opt.method="Nelder-Mead")
summary(m)
summary(m, "DP")
plot(m, which=2, col=4, main="Barolo log10(price)")
# cfr Figure 4.7 of Azzalini & Capitanio (2014), p.107
detach(barolo)
#-----
# examples with multivariate response
#
m3 <- selm(cbind(BMI, LBM) ~ WCC + RCC, family="SN", data=ais)
plot(m3, col=2, which=2)
summary(m3, "dp")
coef(m3)
coef(m3, vector=FALSE)
#
data(wines)
m28 <- selm(cbind(chloride, glycerol, magnesium) ~ 1, family="ST", data=wines)
dp28 <- coef(m28, "DP", vector=FALSE) 
pcp28 <- coef(m28, "pseudo-CP", vector=FALSE) 
\donttest{# the next statement takes a little more time than others
plot(m28)
}
#
m4 <- selm(cbind(alcohol,sugar)~1, family="ST", data=wines)
m5 <- selm(cbind(alcohol,sugar)~1, family="ST", data=wines, fixed=list(alpha=0))
print(1 - pchisq(2*as.numeric(logLik(m4)-logLik(m5)), 2)) # test for symmetry
}

\keyword{regression}
\keyword{univar}
\keyword{multivariate} 
%-------------------------
%% next example has been superseded by introduction of profile.selm
%
% \donttest{
% # example of computation and plot of a (relative twice) profile log-likelihood;
% # since it takes some time, set a coarse grid of nu values
% nu.vector <- seq(3, 8, by=0.5) 
% logL <- numeric(length(nu.vector))
% for(k in 1:length(nu.vector)) { 
%    m28.f <- selm(cbind(chloride, glycerol, magnesium) ~ 1, family="ST", 
%          fixed=list(nu=nu.vector[k]),  data=wines)
%    logL[k] <- logLik(m28.f)
%    cat(format(c(nu.vector[k], logL[k])), "\n")
% }
% plot(nu.vector, 2*(logL-max(logL)), type="b")
% ok <- which.max(logL)
% abline(v=nu.vector[ok], lty=2)
% # compare maximum of this curve with MLE of nu in summary(m28, 'dp')
% }
