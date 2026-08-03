#' @rawNamespace useDynLib(Rsolnp, .registration=TRUE)
#' @keywords internal
#' @importFrom stats rnorm runif
#' @importFrom parallel parLapply clusterExport clusterEvalQ
#' @importFrom numDeriv grad jacobian hessian
#' @importFrom truncnorm rtruncnorm
#' @importFrom Rcpp evalCpp
#' @importFrom utils tail
#' @importFrom stats optim
#' @importFrom future.apply future_lapply
"_PACKAGE"

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
