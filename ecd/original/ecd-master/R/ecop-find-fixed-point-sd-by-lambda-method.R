#' Utility to find the fixed point stdev when lambda is given
#' 
#' This utility finds the fixed point stdev when lambda is given. Other smile related parameters
#' are abstracted away via the closure function \code{fn_get_ld1}.
#' This utility is used primarily to solve the fixed point ATM hypothesis (for VIX option smile).
#' Note that this utility alone is not the full solution. Another utility is needed to match
#' the ATM skew, and the two tails (via mu and epsilon).
#' \code{fn_get_ld1} should have the functional signature of \code{fn_get_ld1(sd, lambda, beta=0)}
#' and returns an ecld object accordingly.
#'
#' @param lambda numeric, the lambda parameter. Must be positive. Default: 3.
#' @param beta numeric, the skewness parameter. Default: 0.
#' @param fn_get_ld1 function, takes stdev, lambda, beta as input, return ld1 object
#'                   via \code{ecop.get_ld_triple}. This closure function encapulates
#'                   mu_plus_ratio, epsilon_ratio, atm_imp_k.
#' @param otype character, option type. Default: "c".
#' @param verbose boolean, print debug message. Default: \code{FALSE}.
#'
#' @return numeric, representing stdev.
#'
#' @keywords fixed-point
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecop.find_fixed_point_sd_by_lambda
#'
### <======================================================================>
"ecop.find_fixed_point_sd_by_lambda" <- function(fn_get_ld1, lambda, beta=0,
                                                 otype="c", verbose=FALSE)
{
    
    get_atm_Q_uniroot <- function(sd) {
        ld1 <- tryCatch(fn_get_ld1(sd, lambda, beta),
                        error = function(cond) { return(NULL) })
        if (is.null(ld1)) return(NaN)
        atm_Q1 <- ecld.fixed_point_atm_Q_right(ld1)
        atm_Q2 <- ecld.fixed_point_atm_Q_left(ld1, otype)
        ecd.mp2f(atm_Q1-atm_Q2)
    }

    ld1_imp <- fn_get_ld1(NaN, lambda, beta) # use atm_imp_k as initial value
    sd_imp <- ld1_imp@mu-ld1_imp@rho
    lower_sd <- 0.25*sd_imp
    upper_sd <- 2*sd_imp
    
    repeat {
        atm_Q_diff <- get_atm_Q_uniroot(lower_sd)
        if (!is.na(atm_Q_diff)) break
        lower_sd <- lower_sd * 1.2
        if (verbose) print(sprintf("Bumped up lower_sd ratio: %.3f", lower_sd/sd_imp))
    }
    
    rs <- uniroot(get_atm_Q_uniroot, lower=lower_sd, upper=upper_sd)
    sd <- rs$root
    ld1 <- fn_get_ld1(sd, lambda, beta)
    Q_fit <- ecld.fixed_point_atm_Q_left(ld1, otype)
    Q_data <- ecld.fixed_point_atm_Q_right(ld1)
    if (verbose) print(sprintf("Q fit: %.3f vs data: %.3f prec: %.4f", Q_fit, Q_data, rs$estim.prec))
    
    return(sd)
}
### <---------------------------------------------------------------------->
