#' The ATM RNO related constants and calculations in fixed point model
#'
#' Computes the small sigma limit of ATM location, rho/stdev,
#' ATM skew of \eqn{Q_c}, and the ratio of lambda to ATM skew 
#' under the RNO measure in the fixed point model.
#'
#' @param lambda numeric the lambda parameter.
#' @param atm_ki numeric optional and experimental, use it as override.
#'               This is for experimental purpose, default is \code{NULL}.
#'               A typical override is the \code{sd/sigma}.
#'
#' @return numeric
#'
#' @keywords ATM
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecld.fixed_point_SN0_atm_ki
#' @export ecld.fixed_point_SN0_atm_ki_sd
#' @export ecld.fixed_point_SN0_rho_sd
#' @export ecld.fixed_point_SN0_skew
#' @export ecld.fixed_point_SN0_lambda_skew_ratio
#'
### <======================================================================>
"ecld.fixed_point_SN0_atm_ki" <- function(lambda) {
    lambda <- ecd.mp2f(lambda) # don't need MPFR. it causes unnecessary complexity. maybe for the future.
    
    if (length(lambda) > 1) {
        object <- ecld(lambda[1], sigma=0.1) # dummy
        rs <- ecld.mclapply(object, lambda, ecld.fixed_point_SN0_atm_ki)
        return(rs)
    }
    
    # process one lambda only
    eq_RN0 <- function(ki) {
        kQ = 1/sqrt(2)
        ki = ki*ecd.mp1 # use MPFR to increase precision
        L1 = sqrt(2)*ki*ecld.ogf_star(ecld(1), kQ)
        L = ecld.ogf_star(ecld(lambda), ki)
        ecd.mp2f(log(L1)-log(L)) # use log to increase precision
    }
    rs <- uniroot(eq_RN0, lower=0.5, upper=20)
    return(rs$root) # ATM ki
}
### <---------------------------------------------------------------------->
#' @rdname ecld.fixed_point_SN0_atm_ki
"ecld.fixed_point_SN0_rho_sd" <- function(lambda) {
    C = sqrt(gamma(3/2*lambda)/gamma(lambda/2))
    -ecld.fixed_point_SN0_atm_ki(lambda)/C
}
### <---------------------------------------------------------------------->
#' @rdname ecld.fixed_point_SN0_atm_ki
"ecld.fixed_point_SN0_atm_ki_sd" <- function() {
    f <- function(x) ecld.fixed_point_SN0_rho_sd(x)+1
    uniroot(f, lower=3.1, upper=3.2)$root # 3.179(1)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.fixed_point_SN0_atm_ki
"ecld.fixed_point_SN0_skew" <- function(lambda, atm_ki=NULL) {
    if (is.null(atm_ki)) atm_ki <- ecld.fixed_point_SN0_atm_ki(lambda)
    C = sqrt(exp(1)*pi)
    A = ecd.erfc(1/sqrt(2))
    B = ecd.gamma(lambda/2, atm_ki^(2/lambda))/gamma(lambda/2)
    skew <- C*(A-B)
    return(skew)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.fixed_point_SN0_atm_ki
"ecld.fixed_point_SN0_lambda_skew_ratio" <- function(lambda, atm_ki=NULL) {
    skew <- ecld.fixed_point_SN0_skew(lambda, atm_ki)
    return(lambda/skew)
}
### <---------------------------------------------------------------------->
