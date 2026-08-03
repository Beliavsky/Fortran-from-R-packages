#' Utility to find the fixed point lambda that matches ATM skew
#' 
#' This utility finds the fixed point lambda from larger lambda to smaller lambda
#' until the calculated ATM skew is smaller than ATM skew from data.
#' It uses \code{ecop.find_fixed_point_sd_by_lambda} to locate stdev.
#' Other smile related parameters are abstracted away via the closure function \code{fn_get_ld1}.
#' This utility is used primarily to solve the fixed point ATM hypothesis (for VIX option smile).
#' Note that this utility alone is not the full solution. Another utility is needed to match
#' the two tails (via mu and epsilon). This utility doesn't handle beta either.
#'
#' @param lambda numeric, the lambda parameter.
#' @param step numeric, increment to decrease lambda.
#' @param fn_get_ld1 function, takes stdev, lambda, beta as input, return ld1 object
#'                   via \code{ecop.get_ld_triple}. This closure function encapulates
#'                   mu_plus_ratio, epsilon_ratio, atm_imp_k.
#' @param atm_skew numeric, ATM skew from data.
#' @param k_atm a vector of numeric, range of log-strike to calculate ATM skew via lm.
#' @param ttm numeric, time to expiration, with 1 representing 1 year (365 days).
#' @param otype character, option type. Default: "c".
#' @param verbose boolean, print debug message. Default: \code{FALSE}.
#' @param msg_prefix character, command line message prefix. Default: "".
#' @param min_lambda numeric, do not try lambda lower than this and return it. Default is 1.1.
#'
#' @return numeric, representing lambda.
#'
#' @keywords fixed-point
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecop.find_fixed_point_lambda_by_atm_skew
#'
### <======================================================================>
"ecop.find_fixed_point_lambda_by_atm_skew" <- function(fn_get_ld1, lambda, step,
                                                atm_skew, k_atm, ttm,
                                                otype="c", verbose=TRUE, msg_prefix="",
                                                min_lambda=1.1)
{
    repeat {
        sd <- ecop.find_fixed_point_sd_by_lambda(fn_get_ld1, lambda, beta=0)
        ld1 <- fn_get_ld1(sd, lambda, beta=0)
        atm_skew_fit <- ecld.op_Q_skew_by_k_lm(ld1, k_atm, otype)/sqrt(2*ttm)
        
        if (verbose) {
            print(paste(msg_prefix,
                        "lambda", lambda,
                        sprintf("sd %.4f", sd),
                        sprintf("atm_skew data: %.3f vs fit: %.3f",
                                atm_skew, atm_skew_fit)
            ))
        }
        
        if (atm_skew_fit <= atm_skew) break
        if (lambda - step < min_lambda) return(lambda)
        lambda <- lambda - step
        if (lambda < 1.1) stop("Unable to find lambda below 1.1")
    }
    return(lambda)
}
    
### <---------------------------------------------------------------------->
