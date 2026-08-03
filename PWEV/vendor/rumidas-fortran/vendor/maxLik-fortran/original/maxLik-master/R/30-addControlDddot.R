
## Method to overwrite parameters of an existing MaxControl object
addControlDddot <- function(x, ...) {
   ## add ... to the control
   dddot <- list(...)
   addControlList(x, dddot)
}

setMethod("maxControl", "MaxControl", addControlDddot)
