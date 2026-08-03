#' The Q operator in option pricing model
#'
#' The Q operator generates the normalized implied volatility \eqn{\sigma_1(k)/\sigma}.
#' \code{cld.op_Q_skew} calculates the skew in Q space by ki and +/- dki/2.
#' \code{cld.op_Q_skew_by_k_lm} calculates the skew in Q space by lm on a vector of k.
#' ki is derived internally from \code{(k-mu-rho)/sigma}.
#' \code{ecld.fixed_point_atm_Q_left} is the left hand side of fixed point ATM hypothesis.
#' \code{ecld.fixed_point_atm_Q_right} is the right hand side of fixed point ATM hypothesis,
#'                                     assuming shift is stored in rho.
#' \code{ecld.fixed_point_atm_ki} is the ATM ki in fixed point ATM hypothesis.
#'                                assuming shift is stored in rho.
#' \code{ecld.fixed_point_shift} is the utility for the standard shift algorithm, -(atm_imp_k - mu).
#'
#' @param object an object of ecld class with built-in \eqn{\rho, \epsilon}
#' @param ki numeric, a vector of \eqn{\sigma}-normalized log-strike
#' @param dki numeric, delta of ki for calculating slope
#' @param k numeric, a vector of log-strike
#' @param atm_imp_k numeric, the ATM implied log-strike. It is derived from ATM volatility
#'                  times sqare root of time to expiration.
#' @param otype character, specifying option type:
#'              \code{c} (default) or \code{p}.
#'
#' @return a numeric vector, representing Q or skew of Q.
#'         For \code{ecld.fixed_point_atm_ki}, it is ATM ki.
#'         For \code{ecld.fixed_point_shift}, it is the shift.
#'
#' @keywords Q
#'
#' @author Stephen H. Lihn
#'
#' @export ecld.op_Q
#' @export ecld.op_Q_skew
#' @export ecld.op_Q_skew_by_k_lm
#' @export ecld.fixed_point_atm_Q_left
#' @export ecld.fixed_point_atm_Q_right
#' @export ecld.fixed_point_atm_ki
#' @export ecld.fixed_point_shift
#'
### <======================================================================>
"ecld.op_Q" <- function(object, ki, otype="c")
{
    ecld.validate(object, sged.allowed=TRUE)
    if (is.na(object@mu_D)) object@mu_D <- ecld.mu_D(object)

    k <- ki*object@sigma + object@mu
    L0 <- ecld.ogf(object, k, otype=otype, RN=FALSE)
    epsilon <- if(is.na(object@epsilon)) 0 else object@epsilon
    L <- L0 + epsilon
    ecd.mp2f(ecld.op_V(L, k, otype=otype)/object@sigma)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.op_Q
"ecld.op_Q_skew" <- function(object, ki, dki=0.1, otype="c") {
    ecld.validate(object, sged.allowed=TRUE)

    if (length(ki) > 1) {
        S <- ecld.sapply(object, ki, function(x) ecld.op_Q_skew(object, x, dki, otype))
        return(ecd.mp2f(S))
    }
    Q <- ecld.op_Q(object, c(ki-dki/2, ki+dki/2), otype)
    ecd.mp2f((Q[2]-Q[1])/dki) # return a negative slope
}
### <---------------------------------------------------------------------->
#' @rdname ecld.op_Q
"ecld.op_Q_skew_by_k_lm" <- function(object, k, otype="c") {
    ecld.validate(object, sged.allowed=TRUE)

    ki <- (k - object@mu - object@rho)/object@sigma
    Q <- ecld.op_Q(object, ki, otype)
    Q.lm = lm(Q ~ ki)
    Q.c = Q.lm$coefficients
    unname(Q.c[2]) # slope
}
### <---------------------------------------------------------------------->
#' @rdname ecld.op_Q
"ecld.fixed_point_atm_Q_left" <- function(object, otype="c") {
    ecld.validate(object, sged.allowed=FALSE)
    
    atm_ki <- ecld.fixed_point_atm_ki(object)
    ecld.op_Q(object, atm_ki, otype)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.op_Q
"ecld.fixed_point_atm_ki" <- function(object) {
    ecld.validate(object, sged.allowed=FALSE)

    atm_imp_k <- object@mu - object@rho
    ecd.mp2f((atm_imp_k - 2*object@mu)/object@sigma)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.op_Q
"ecld.fixed_point_atm_Q_right" <- function(object) {
    ecld.validate(object, sged.allowed=FALSE)
    
    atm_imp_k <- object@mu - object@rho
    ecd.mp2f(sqrt(2)*atm_imp_k/object@sigma)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.op_Q
"ecld.fixed_point_shift" <- function(object, atm_imp_k) {
    -ecd.mp2f(atm_imp_k - object@mu)
}
### <---------------------------------------------------------------------->

