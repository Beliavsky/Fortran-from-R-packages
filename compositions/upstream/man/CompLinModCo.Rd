\name{CompLinModCoReg}
\alias{CompLinModCoReg}
\title{Compositional Linear Model of Coregionalisation}
\description{
  Creates a Variogram model according to the linear model of spatial
  corregionalisation for a compositional geostatistical analysis.
}
\usage{
CompLinModCoReg(formula,comp,D=ncol(comp),envir=environment(formula))
          }
\arguments{
  \item{formula}{A formula without left side providing a formal
    description of a variogram model.}
  \item{comp}{a compositional dataset, needed to provide the frame size}
  \item{D}{The dimension of the multivariate dataset}
  \item{envir}{The enviroment in which formula should be interpreted.}
}
\value{
  A variogram function, with the extra class "\code{CompLinModCoReg}". 
}
\details{
  The linear model of coregionalisation uses the fact that sums of
  valid variogram models are valid variograms, and that scalar variograms
  multiplied with a positive definite matrix are valid variograms for
  vector-valued random functions. 
  
  This command computes such a variogram function from a
  formal description, via a formula without left-hand side. 
  The right-hand side of the formula is a sum. Each summand is either a
  product of a matrix description and a scalar variogram description or
  only a scalar variogram description. Scalar variogram descriptions are
  either formal function calls to
  \describe{
    \item{\code{sph(range)}}{ for spherical variogram with range \code{range}}
    \item{\code{exp(range)}}{ for an exponential variogram with effective range 
      \code{range}}
    \item{\code{gauss(range)}}{ for a Gaussian variogram with effective range
    \code{range}}
    \item{\code{gauss(range)}}{ for a cardinal sine variogram with range parameter
    \code{range}}
    \item{\code{pow(range)}}{ for an power variogram with range parameter
      \code{range}}
    \item{\code{lin(unit)}}{ linear variogram 1 at \code{unit}.}
    \item{\code{nugget()}}{ for adding a nuggeteffect.}
  }
  Alternatively it can be any expression, which will be evaluated in
  envir and should depende on a dataset of distantce vectrs \code{h}.
  An effective range is that distance at which one reaches the sill (for spherical)
  of 95\% of its values (for all other models). Parametric ranges are given for those
  models that do not have an effective range formula.
  \cr
  The matrix description always comes first. It can be \code{R1} for a
  rank 1 matrix; \code{PSD} for a Positive Semidefinite matrix; \(S\)
  for a scalar Sill factor to be multiplied with the identity matrix; or any other
  construct evaluating to a matrix, like e.g. a function of some parameters with
  default values, that if called is evaluated to a positive semidefinite
  matrix. \code{R1} and \code{PSD} can also be written as calls
  providing a vector or respectively a matrix providing the parameter.
  \cr
  The variogram is created with default parameter values. The parameters
  can later be modified by modifiying the default parameter with
  assignments like \code{formals(vg)$sPSD1 =
    parameterPosdefMat(4*diag(5))}.
  We would anyway expect you to fit the model to the data by a command
  like \code{fit.lmc(logratioVariogram(...),CompLinModCoReg(...))}  
}
\references{
  What to cite??
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{vgram2lrvgram}},
  \code{\link{vgmFit2lrv}}
}
\examples{
\dontrun{
data(juraset)
X <- with(juraset,cbind(X,Y))
comp <- acomp(juraset,c("Cd","Cu","Pb","Co","Cr"))
CompLinModCoReg(~nugget()+sph(0.5)+R1*exp(0.7),comp)
CompLinModCoReg(~nugget()+R1*sph(0.5)+R1*exp(0.7)+(0.3*diag(5))*gauss(0.3),comp)
CompLinModCoReg(~nugget()+R1*sph(0.5)+R1(c(1,2,3,4,5))*exp(0.7),comp)
}
}
\keyword{multivariate}
