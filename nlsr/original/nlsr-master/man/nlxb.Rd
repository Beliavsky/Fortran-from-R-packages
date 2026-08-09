\name{nlxb}
\alias{nlxb}
\title{nlxb: nonlinear least squares modeling by formula}
\usage{
nlxb(
  formula,
  data = parent.frame(),
  start,
  trace = FALSE,
  lower = NULL,
  upper = NULL,
  weights = NULL,
  control = list(),
  ...
)
}
\arguments{
\item{formula}{The modeling formula. Looks like 'y~b1/(1+b2*exp(-b3*T))'}

\item{data}{a data frame containing data for variables
used in the formula that are NOT the parameters. This
data may also be defined in the parent frame i.e.,
'global' to this function}

\item{start}{MUST be a named vector with all elements present
e.g., start=c(b1=200, b2=50, b3=0.3)}

\item{trace}{TRUE for console output during execution}

\item{lower}{a vector of lower bounds on the parameters. 
If a single number, this will be applied to all parameters
Default \code{NULL}.}

\item{upper}{a vector of upper bounds on the parameters. If a single number, 
this will be applied to all parameters. Default \code{NULL}.}

\item{weights}{A vector of fixed weights or a function or formula
producing one.  See the Details below.}

\item{control}{a list of control parameters. See nlsr.control().}

\item{...}{additional data needed to evaluate the modeling functions}
}
\value{
list of solution elements
  \item{resid}{weighted residuals at the proposed solution}
  \item{jacobian}{Jacobian matrix at the proposed solution}
  \item{feval   }{residual function evaluations used to reach solution from starting parameters}
  \item{jeval   }{Jacobian function (or approximation) evaluations used}
  \item{coefficients}{a named vector of proposed solution parameters}
  \item{ssquares}{weighted sum of squared residuals (often the deviance)}
  \item{lower   }{lower bounds on parameters}
  \item{upper   }{upper bounds on parameters}
  \item{maskidx }{vector if indices of fixed (masked) parameters}
  \item{weights0}{weights specified in function call}
  \item{weights }{weights at the final solution}
  \item{formula }{the modeling formula}
  \item{resfn   }{the residual function (unweighted) based on the formula}
}
\description{
A simplified and hopefully robust alternative to finding
the nonlinear least squares minimizer that causes
'formula' to give a minimal residual sum of squares.
}
\details{
nlxb is particularly intended to allow for the
resolution of very ill-conditioned or else near
zero-residual problems for which the regular nls()
function is ill-suited. 

This variant uses a qr solution without forming the sum
of squares and cross products t(J)%*%J

Neither this function nor \code{nlfb} have provision for parameter
scaling (as in the \code{parscale} control of \code{optim} and
package \code{optimx}). This would be more tedious than difficult to
introduce, but does not seem to be a priority feature to add.


  There are many controls, and some of them are important for \code{nlxb}.
  In particular, if the derivatives needed for developing the Jacobian are
  NOT in the derivatives table, then we must supply code elsewhere as 
  specified by the control \code{japprox}. This was originally just for
  numerical approximations, with the character strings "jafwd", "jaback",
  "jacentral" and "jand" leading to the use of a forward, backward, central
  or package \code{numDeriv} approximation. However, it is also possible to
  use code embedded in the residual function created using the \code{formula}.
  This is particularly useful for \code{selfStart} models, and we use the
  character string "SSJac" to point to such Jacobian code. Note how the
  starting parameter vector is found using the \code{getInitial} function
  from the \code{stats} package as in an example below.
  
The \code{weights} argument can be a vector of fixed weights, in which
case the objective function that will be 
minimized is the sum of squares where each residual is multiplied by the 
square root of the corresponding weight. Default \code{NULL} implies 
unit weights.  

\code{weights} may alternatively be a function with header \code{function(parms, resids)} to compute such a vector, or a formula
whose right hand side gives an expression for the weights.  Variables
in the expression may include the following:
\describe{
\item{A variable named \code{resid}}{The current residuals.}
\item{A variable named \code{fitted}}{The right hand side of the model formula.}
\item{Parameters}{The parameters of the model.}
\item{Data}{Values from \code{data}.}
\item{Vars}{Variables in the environment of the formula.}
}
}
\examples{
library(nlsr)
weed <- c(5.308, 7.24, 9.638, 12.866, 17.069, 23.192, 31.443,
        38.558, 50.156, 62.948, 75.995, 91.972)
tt <- 1:12
weeddf <- data.frame(tt, weed)
frm <- 
wmodu <- weed ~ b1/(1+b2*exp(-b3*tt)) # Unscaled
## nls from unit start FAILS
start1<-c(b1=1, b2=1, b3=1)
hunls1 <- try(nls(wmodu, data=weeddf, start=start1, trace=FALSE))
if (! inherits(hunls1, "try-error")) print(hunls1) ## else cat("Failure -- try-error\n")
## nlxb from unit start
hunlx1 <- try(nlxb(wmodu, data=weeddf, start=c(b1=1, b2=1, b3=1), trace=FALSE))
if (! inherits(hunlx1, "try-error")) print(hunlx1)

st2h<-c(b1=185, b2=10, b3=.3)
#' hunls2 <- try(nls(wmodu, data=weeddf, start=st2h, trace=FALSE))
if (! inherits(hunls1, "try-error")) print(hunls1) ## else cat("Failure -- try-error\n")
## nlxb from unit start
hunlx1 <- try(nlxb(wmodu, data=weeddf, start=st2h, trace=FALSE))
if (! inherits(hunlx1, "try-error")) print(hunlx1)

# Functional models need to use a Jacobian approximation or external calculation.
# For example, the SSlogis() selfStart model from \code{stats} package.

# nls() needs NO starting value
hSSnls<-try(nls(weed~SSlogis(tt, Asym, xmid, scal), data=weeddf))
summary(hSSnls)
# We need to get the start for nlxb explicitly
stSS <- getInitial(weed ~ SSlogis(tt, Asym, xmid, scal), data=weeddf)
hSSnlx<-try(nlxb(weed~SSlogis(tt, Asym, xmid, scal), data=weeddf, start=stSS))
hSSnlx

# nls() can only bound parameters with algorithm="port"
#   and minpack.lm is unreliable in imposing bounds, but nlsr copes fine.
lo<-c(0, 0, 0)
up<-c(190, 10, 2) # Note: start must be admissible.
bnls0<-try(nls(wmodu, data=weeddf, start=st2h,
         lower=lo, upper=up)) # should complain and fail
 
bnls<-try(nls(wmodu, data=weeddf, start=st2h,
         lower=lo, upper=up, algorith="port"))
summary(bnls)
bnlx<-try(nlxb(wmodu, data=weeddf, start=st2h, lower=lo, upper=up))
bnlx

# nlxb() can also MASK (fix) parameters. The mechanism of maskidx from nls
# is NO LONGER USED. Instead we set upper and lower parameters equal for
# the masked parameters. The start value MUST be equal to this fixed value.
lo<-c(190, 0, 0) # mask first parameter
up<-c(190, 10, 2)
strt <- c(b1=190, b2=1, b3=1)
mnlx<-try(nlxb(wmodu, start=strt, data=weeddf, 
         lower=lo, upper=up))
mnlx
mnls<-try(nls(wmodu, data=weeddf, start=strt,
         lower=lo, upper=up, algorith="port"))
summary(mnls)

# Try first parameter masked and see if we get SEs 
lo<-c(200, 0, 0) # mask first parameter
up<-c(100, 10, 5)
strt <- c(b1=200, b2=1, b3=1)
mnlx<-try(nlxb(wmodu, start=strt, data=weeddf, 
         lower=lo, upper=up))
mnlx
mnls<-try(nls(wmodu, data=weeddf, start=strt,
         lower=lo, upper=up, algorith="port"))
summary(mnls) 

# Try with weights on the observations
mnlx<-try(nlxb(wmodu, start=strt, data=weeddf, 
               weights = ~ 1/weed))
mnlx

}
\author{
J C Nash 2014-7-16   nashjc _at_ uottawa.ca
}
