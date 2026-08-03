#' Get triple list of ecld objects by stdev
#' 
#' Construct triple list of ecld objects by stdev, with lambda, and ratios related to stdev.
#' This utility is used primarily in fixed point ATM hypothesis (when simulating VIX option smile).
#'
#' @param lambda numeric, the lambda parameter. Must be positive. Default: 3.
#' @param sd numeric, the stdev parameter. Must be positive. Default: 1.
#' @param beta  numeric, the skewness parameter. Default: 0.
#' @param mu_plus_ratio numeric, numeric, excess value in addition to \code{mu_D},
#'                      relative to the stdev. Default: 0.
#' @param epsilon_ratio numeric, epsilon ratio relative to the stdev. Default: 0.
#' @param atm_imp_k numeric, ATM implied log-strike. It is derived from ATM volatility
#'                  times sqare root of time to expiration. If provided, it is used to calculate
#'                  the fixed point shift, -(atm_imp_k - mu). Default: NaN.
#'                  the rho slot in ld1 is populated with the value.
#' @param fn_shift function, takes an ecld object and return the fixed point shift, -(atm_imp_k - mu).
#'                 the rho slot in ld1 is populated with the value from this function.
#'                 This serves as secondary method if you don't want to provide atm_imp_k directly.
#'
#' @return a triple list of ecld objects.
#' ld0 has mu=0 as vanila object; ld1 has mu and rho as prescribed; ld2 has mu=mu_D.
#'
#' @keywords constructor
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecop.get_ld_triple
#'
#' @examples
#' lds <- ecop.get_ld_triple(3, 0.1)
#' ld1 <- lds$ld1
### <======================================================================>
"ecop.get_ld_triple" <- function(lambda = 3, sd = 1, beta = 0,
                    mu_plus_ratio = 0, epsilon_ratio = 0,
                    atm_imp_k = NaN, fn_shift = NULL)
{
    if (is.na(sd) & !is.na(atm_imp_k)) sd <- atm_imp_k
    if (is.na(sd)) stop("Error: sd is not provided")

    ld0 <- ecld.from_sd(lambda, sd, beta=beta)
    sigma <- ld0@sigma

    sd0 <- ecld.sd(ld0)
    mu_plus <- mu_plus_ratio * sd0
    mu_D <- ecld.mu_D(ld0)
    ld0@mu_D <- mu_D

    ld1 <- ecld(lambda=lambda, sigma=sigma, beta=beta, mu=mu_D+mu_plus)
    ld1@mu_D <- mu_D
    ld1@epsilon <- epsilon_ratio * sd0
    if (! is.na(atm_imp_k)) ld1@rho <- ecld.fixed_point_shift(ld1, atm_imp_k)
    if (! is.null(fn_shift))  ld1@rho <- fn_shift(ld1)

    ld2 <- ecld(lambda=lambda, sigma=sigma, beta=beta, mu=mu_D)
    ld2@mu_D <- mu_D
    ld2@epsilon <- epsilon_ratio * sd0
    
    list(ld0=ld0, ld1=ld1, ld2=ld2)
}
### <---------------------------------------------------------------------->
