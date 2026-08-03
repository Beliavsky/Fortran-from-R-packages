### Methods for accessing loglik value maximum likelihood estimates

logLik.summary.maxLik <- function( object, ...) {
   ll <- object$loglik
   attr(ll, "df") <- sum(activePar(object))
   ll
}

logLik.maxLik <- function( object, ...) {
   ll <- maxValue(object)
   attr(ll, "df") <- sum(activePar(object))
   ll
}
