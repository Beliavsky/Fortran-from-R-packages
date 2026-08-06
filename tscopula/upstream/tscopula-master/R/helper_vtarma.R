#' Conditional density of VT-ARMA process
#'
#' @param x vector of points at which density should be calculated.
#' @param tscmmod an object of class \linkS4class{tscmfit} based on underlying copula
#' of class \linkS4class{armacopula}.
#' @param armamean conditional mean of underlying ARMA process.
#' @param armasd conditional standard deviation of underlying ARMA process.
#'
#' @return vector with same length as x.
#' @keywords internal
#'
dcondvtarma <- function(x, tscmmod, armamean, armasd){
  margmod <- tscmmod@margin
  u <- pmarg(margmod, x)
  v <- vtrans(tscmmod@tscopula@Vtransform, u)
  dv <- dnorm(qnorm(v),
              mean = armamean,
              sd = armasd,
              log = TRUE) - dnorm(qnorm(v), log = TRUE) + dmarg(margmod, x, log = TRUE)
  dv[is.nan(dv)] <- -Inf
  dv[is.nan(dv)] <- Inf
  return(exp(dv))}

#' Conditional distribution function of VT-ARMA Process
#'
#' @param q point at which CDF should be calculated.
#' @param tscmmod an object of class \linkS4class{tscmfit} based on underlying copula
#' of class \linkS4class{armacopula}.
#' @param armamean conditional mean of underlying ARMA process.
#' @param armasd conditional standard deviation of underlying ARMA process.
#'
#' @return a scalar value.
#' @keywords internal
#'
pcondvtarma <- function(q, tscmmod, armamean, armasd) {
  integrate(
    function(t)
      dcondvtarma(t, tscmmod, armamean, armasd),
    lower = -Inf,
    upper = q
  )$value
}

#' Conditional quantiles of VT-ARMA process
#'
#' @param p point at which quantile should be calculated.
#' @param tscmmod an object of class \linkS4class{tscmfit} based on underlying copula
#' of class \linkS4class{armacopula}.
#' @param armamean conditional mean of underlying ARMA process.
#' @param armasd conditional standard deviation of underlying ARMA process.
#'
#' @return a scalar value.
#' @keywords internal
#'
qcondvtarma <- function(p, tscmmod, armamean, armasd) {
  uniroot(function(t)
    pcondvtarma(t, tscmmod, armamean, armasd) - p,
    lower = -30,
    upper = 30)$root
}

pcondvtarma <- Vectorize(pcondvtarma, c("q", "armamean", "armasd"))
qcondvtarma <- Vectorize(qcondvtarma, c("p", "armamean", "armasd"))

#' Quantile calculation method for VT-ARMA models
#'
#' @param x an object of class \linkS4class{tscmfit} based on underlying copula
#' of class \linkS4class{armacopula}.
#' @param alpha a scalar probability value
#' @param last logical value asserting that only the last volatility
#' prediction should be returned
#'
#' @return a vector of the same length as the data embedded in the tscmfit object.
#' @export
#'
setMethod("quantile", c(x = "tscmfit"), function(x, alpha, last = FALSE){
  margmod <- x@margin
  tscopula <- x@tscopula
  U <- pmarg(margmod, x@data)
  if (!(is(tscopula, "vtscopula")))
    stop("tscopula must be vtscopula")
  vt <- tscopula@Vtransform
  V <- vtrans(vt, U)
  armacop <- tscopula@Vcopula
  if (!(is(armacop, "armacopula")))
    stop("Underlying copula must be ARMA")
  series <- kfilter(armacop, V)
  mu_t <- series[, "mu_t"]
  sigma_t <- series[, "sigma_t"]
  if (last){
    mu_t <- mu_t[length(V)]
    sigma_t <- sigma_t[length(V)]
  }
  VaR <- qcondvtarma(alpha,
                     tscmmod = x,
                     armamean = mu_t,
                     armasd = sigma_t)
  if (!last)
    attributes(VaR) <- attributes(x@data)
  VaR
}
)
