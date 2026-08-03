#' The ATM volatility and skew of \eqn{Q_p} in quartic model
#'
#' Compute the ATM location and ATM skew of \eqn{Q_p} in quartic model.
#'
#' @param object an object of ecd class
#' @param ki numeric, order of the moment to be computed
#' @param atm_ki numeric, if provided, take it as is without calculating again
#' @param dki numeric, delta of ki for calculating slope
#' @param otype character, specifying option type with either c or p.
#' @param lower numeric, optional value to specify the lower bound of ATM root finding.
#'              This is often needed when the smile is collapsed in the left wing.
#' @param upper numeric, optional value to specify the upper bound of ATM root finding
#'              This is often needed when the smile is collapsed significantly in the right wing.
#'
#' @return numeric
#'
#' @keywords ATM
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecld.quartic_Qp
#' @export ecld.quartic_Q
#' @export ecld.quartic_Qp_atm_ki
#' @export ecld.quartic_Qp_rho
#' @export ecld.quartic_Qp_atm_skew
#' @export ecld.quartic_Qp_skew
#'
#' @examples
#' \dontrun{
#'     ld <- ecld.quartic(sigma=0.001*ecd.mp1, epsilon=0, rho=0, mu_plus=0)
#'     ecld.quartic_Qp_atm_ki(ld, lower=-12, upper=-11)
#'     ecld.quartic_Qp_atm_skew(ld, lower=-12, upper=-11)
#' }
### <======================================================================>
"ecld.quartic_Qp" <- function(object, ki) {
    ecld.quartic_Q(object, ki, otype="p")
}
### <---------------------------------------------------------------------->
#' @rdname ecld.quartic_Qp
"ecld.quartic_Q" <- function(object, ki, otype) {
    rho <- if(is.na(object@rho)) 0 else object@rho
    k = ki*object@sigma + object@mu + rho
    IV <- ecd.mp2f(ecld.op_VL_quartic(object, k, otype=otype))
    ecd.mp2f(IV/object@sigma)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.quartic_Qp
"ecld.quartic_Qp_atm_ki" <- function(object, lower=-50, upper=-1.37)
{
    target <- sqrt(240)
    f <- function(ki) ecld.quartic_Qp(object, ki)-target
    while (lower > -200 & is.na(f(lower))) {
        upper = lower # this will definitely mess up in a collapsed call smile
        lower = lower*1.5
    }
    while (lower > -200 & f(lower) < 0) {
        upper = lower
        lower = lower*1.5
    }
    while (upper > -200 & is.na(f(upper)) ) upper = upper*1.5 # bump it up
    while (upper <= -1 & f(upper) > 0 ) upper = upper/1.5 # bump it down
    
    if (f(lower) * f(upper) < 0 & upper > lower) {
        rs <- uniroot(f, lower=lower, upper=upper)
        return(rs$root)
    }
    return(NaN)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.quartic_Qp
"ecld.quartic_Qp_rho" <- function(object, atm_ki=NaN, lower=-50, upper=-1.37)
{
    if (is.na(atm_ki)) {
      atm_ki <- ecld.quartic_Qp_atm_ki(object, lower, upper)
    }
    -ecd.mp2f(object@mu + atm_ki * object@sigma)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.quartic_Qp
"ecld.quartic_Qp_skew" <- function(object, ki, dki=0.1) {
    if (length(ki) > 1) {
       S <- ecld.sapply(object, ki, function(x) ecld.quartic_Qp_skew(object, x, dki))
       return(ecd.mp2f(S))
    }
    Q <- ecld.quartic_Qp(object, c(ki, ki+dki))
    ecd.mp2f((Q[2]-Q[1])/dki) # return a negative slope
}
### <---------------------------------------------------------------------->
#' @rdname ecld.quartic_Qp
"ecld.quartic_Qp_atm_skew" <- function(object, dki=0.1, lower=-50, upper=-1.37) {
    atm_ki <- ecld.quartic_Qp_atm_ki(object, lower, upper)
    ecld.quartic_Qp_skew(object, atm_ki, dki)
}
### <---------------------------------------------------------------------->
