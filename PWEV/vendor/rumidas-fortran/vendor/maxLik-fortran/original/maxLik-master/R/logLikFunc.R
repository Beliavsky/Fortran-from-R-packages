if( getRversion() >= "2.15.1" ) {
   globalVariables( c( "lastFuncGrad", "lastFuncParam" ) )
}

## objective function:
## sum over possible individual likelihoods
logLikFunc <- function(theta,
                       fnOrig,  # the original user-supplied function we wrap here
                       gradOrig,
                       hessOrig,
                           # Arguments "gradOrig" and "hessOrig" are just for compatibility with
                           # logLikGrad() and logLikHess()
                       start = NULL, fixed = NULL, sumObs = TRUE, ...) {
   if(missing(fnOrig)) {
      stop("Cannot compute the objective function value: no objective function supplied")
   }
   theta <- addFixedPar( theta = theta, start = start, fixed = fixed, ...)
   result <- fnOrig( theta, ... )

   ## save gradients and the corresponding parameter values
   assign( "lastFuncGrad", attr( result, "gradient" ), inherits = TRUE )
   assign( "lastFuncParam", theta, inherits = TRUE )

   if( sumObs ) {
      result <- sumKeepAttr( result )
      g <- attributes( result )$gradient
      if( !is.null( g ) ) {
         g <- sumGradients( g, length( theta ) )
         names( g ) <- names( theta )
         if( !is.null( fixed ) ) {
            g <- g[ !fixed ]
         }
         attributes( result )$gradient <- g
      }
   }

   return( result )
}
