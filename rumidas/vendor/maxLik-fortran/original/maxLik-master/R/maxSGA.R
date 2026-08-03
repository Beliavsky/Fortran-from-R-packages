
maxSGA <- function(fn=NULL, grad=NULL, hess=NULL, start,
                   nObs,
                   constraints=NULL,
                   finalHessian=FALSE,
                   fixed=NULL,
                   control=NULL,
                   ...) {
   ## Newton-Raphson maximisation
   ## Parameters:
   ## fn          - the function to be maximized.  Returns either scalar or
   ##               vector value with possible attributes 
   ##               constPar and newVal
   ## grad        - gradient function (numeric used if missing).  Must return either
   ##               * vector, length=nParam
   ##               * matrix, dim=c(nObs, 1).  Treated as vector
   ##               * matrix, dim=c(M, nParam), where M is arbitrary.  In this case the
   ##                 rows are simply summed (useful for maxBHHH).
   ## hess        - hessian function (used only for finalHessian, otherwise ignored)
   ## start       - initial parameter vector (eventually w/names)
   ## ...         - extra arguments for fn()
   ## finalHessian  include final Hessian?  As computing final hessian does not carry any extra penalty for NR method, this option is
   ##               mostly for compatibility reasons with other maxXXX functions.
   ##               TRUE/something else  include
   ##               FALSE                do not include
   ## fixed         index vector, which parameters to keep fixed
   ##
   ## RESULTS:
   ## a list of class "maxim":
   ## maximum     function value at maximum
   ## estimate    the parameter value at maximum
   ## gradient        gradient
   ## hessian         Hessian
   ## code        integer code of success:
   ##             1 - gradient close to zero
   ##             2 - successive values within tolerance limit
   ##             3 - could not find a higher point (step error)
   ##             4 - iteration limit exceeded
   ##             100 - initial value out of range
   ## message     character message describing the code
   ## iterations  number of iterations
   ## type        "Newton-Raphson maximisation"
   ##
   ## ------------------------------
   ## Add parameters from ... to control
   if(!inherits(control, "MaxControl")) {
      mControl <- addControlList(maxControl(gradtol=0), control)
   }
   else {
      mControl <- control
   }
   mControl <- addControlList(mControl, list(...), check=FALSE)
   ##
   argNames <- c(c("fn", "grad", "hess", "start",
                   "fixed", "control"),
                 openParam(mControl))
                           # Here we allow to submit all parameters outside of the
                           # 'control' list.  May eventually include only a
                           # subset here
   ## ensure that 'fn', 'grad', and 'hess' do not take any arguments that maxSGA eats up
   if(!is.null(fn)) {
      checkFuncArgs( fn, argNames, "fn", "maxSGA" )
   }
   if( !is.null( grad ) ) {
      checkFuncArgs( grad, argNames, "grad", "maxSGA" )
   }
   if( !is.null( hess ) ) {
      checkFuncArgs( hess, argNames, "hess", "maxSGA" )
   }
   ## ensure that at least 'fn' or 'grad' are supplied
   if(is.null(fn) & is.null(grad)) {
      stop("maxSGA/maxAdam requires at least 'fn' or 'grad' to be supplied")
   }
   if(length(start) < 1) {
      stop("'start' must be of positive length!")
   }
   ## establish the active parameters.  Internally, we just use 'activePar'
   fixed <- prepareFixed( start = start, activePar = NULL,
      fixed = fixed )
   ## chop off the control args from ... and forward the new ...
   dddot <- list(...)
   dddot <- dddot[!(names(dddot) %in% openParam(mControl))]
   cl <- list(start=start,
              finalHessian=finalHessian,
              fixed=fixed,
              control=mControl,
              optimizer="SGA")
   if(length(dddot) > 0) {
      cl <- c(cl, dddot)
   }
   ##
   if(is.null(constraints)) {
      ## call maxSGACompute with the modified ... list
      cl <- c(quote(maxSGACompute),
              fn=logLikFunc, grad=logLikGrad, hess=logLikHess,
              fnOrig = fn, gradOrig = grad, hessOrig = hess,  # these are forwarded to the logLikAttr
              nObs=nObs,
              cl)
      result <- eval(as.call(cl))
   } else {
      if(identical(names(constraints), c("ineqA", "ineqB"))) {
         stop("Inequality constraints not implemented for maxSGA")
      } else if(identical(names(constraints), c("eqA", "eqB"))) {
                           # equality constraints: A %*% beta + B = 0
         cl <- c(quote(sumt),
                 fn=fn, grad=grad, hess=hess,
                 maxRoutine=maxSGA,
                 constraints=list(constraints),
                 cl)
         result <- eval(as.call(cl))
      } else {
         stop("maxNR only supports the following constraints:\n",
              "constraints=list(ineqA, ineqB)\n",
              "\tfor A %*% beta + B >= 0 linear inequality constraints\n",
              "current constraints:",
              paste(names(constraints), collapse=" "))
      }
   }
   ## Save the objective function
   result$objectiveFn <- fn
   ##
   return( result )
}
