\name{vgmFit}
\alias{vgmFit}
\alias{vgmFit2lrv}
\alias{vgmGetParameters}
\alias{vgmSetParameters}
\alias{vgmGof}
\alias{fit.lmc}
\alias{fit.lmc.logratioVariogram}
\title{Compositional variogram model fitting}
\description{
  Fits a parametric variogram model to an empirical logratio-Variogram 
}
\usage{
vgmFit2lrv(emp,vg,...,mode="log",psgn=rep(-1,length(param)),print.level=1)
\method{fit.lmc}{logratioVariogram}(v,model,...,mode="log",psgn=rep(-1,length(param)),print.level=1)
vgmFit(emp,vg,...,mode="log",psgn=rep(-1,length(param)),print.level=1)
vgmGof(p = vgmGetParameters(vg), emp, vg, mode = "log")
vgmGetParameters(vg,envir=environment(vg))
vgmSetParameters(vg,p)
fit.lmc(v,...)
          }
\arguments{
  \item{emp}{An empirical logratio-Variogram as e.g. returned by \code{\link{logratioVariogram}}}
  \item{v}{An empirical logratio-Variogram as e.g. returned by \code{\link{logratioVariogram}}}
  \item{vg}{A compositional clr-variogram (or ilt-vagriogram) model function.}
  \item{model}{A compositional clr-variogram (or ilt-vagriogram) model function, output of a call to
  .}
  \item{\dots}{further parameters to \code{\link{nlm}}}
  \item{mode}{either "ls" or "log"  for selection of either using either
    least squares or least squares on logarithmic values.}
  \item{psgn}{Contains a parameter code for each of the parameters. -1
    means the parameter should be used as is. 0 means the parameter is
    nonnegativ and 1 means the parameter is striktly positiv. This
    allows to provide parameter limits if the fitting procedure fails.}
  \item{print.level}{The print.level of \code{\link{nlm}}. 0 for no
    printing. 1 for a rough information about the sucess and 2 for step
    by step printing.}
  \item{p}{Is the parameter of the variogram model in linearized form as
    e.g. 
    returned by \code{vgmGetParameters}.
  }
  \item{envir}{The environment the default parameters of the model
    should be evaluated in.}
}
\value{
  
  \code{vgmFit2lrv} returns a list of two elements.
    \item{nlm}{The result of \code{\link{nlm}} containing covergence
      codes.}
    \item{vg}{A version of \code{vg} but with default parameters
      modified according to the fitting.}
  \code{vgmGof} returns a scalar quantifiying the goodness of fit, of a
  model and an empirical variogram.
  \cr
  \code{vgmGetParameters} extracts the default values of a variogram model
  function to a parameter vector. It returns a numeric vector.
  \cr
  \code{vgmSetParameters} does the inverse operation and modifies the
  default according to the new values in \code{p}. It returns \code{vg}
  with modifiend default parameter values.
}
\details{
  The function is mainly a wrapper to \code{\link{nlm}} specifying the
  an objective function for modell fitting, taking the starting values
  of fitting procedure from the default arguments and writing the
  results back. Variogram model fitting is more an art than a straight
  forward procedure. Fitting procedures typically only find a right
  optimum if reasonable starting parameters are provided.
  The fit should
  be visually checked afterwards.
  \cr
  The meaning of \code{psgn} is subject to change. We will probably
  provide a more automatic procedure later.
  \cr
  \code{vgmFit} is a copy of \code{vgmFit2lrv}, but deprecated. The name
  will later be used for other functionality. 
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{vgram2lrvgram}},
  \code{\link{CompLinModCoReg}},
  \code{\link{logratioVariogram}}
}
\examples{
\dontrun{
data(juraset)
X <- with(juraset,cbind(X,Y))
comp <- acomp(juraset,c("Cd","Cu","Pb","Co","Cr"))
lrv <- logratioVariogram(comp,X,maxdist=1,nbins=10)
fff <- CompLinModCoReg(~nugget()+sph(0.5)+R1*exp(0.7),comp)
fit <- vgmFit(lrv,fff)
fit
fff(1:3)
plot(lrv,lrvg=vgram2lrvgram(fit$vg))
}
}
\keyword{multivariate}
