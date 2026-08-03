#' The ATM RNO related constants and calculations in quartic model
#'
#' Computes the small sigma limit of ATM location, rho/stdev,
#' and ATM skew of \eqn{Q_p} under the RNO measure in quartic model.
#' Computes the maximum risk-neutral violation as an extension of RN0 measure.
#'
#' @param sigma numeric, the volatility parameter
#'
#' @return numeric
#'
#' @keywords ATM
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecld.quartic_SN0_atm_ki
#' @export ecld.quartic_SN0_rho_stdev
#' @export ecld.quartic_SN0_skew
#' @export ecld.quartic_SN0_max_RNV
#'
### <======================================================================>
"ecld.quartic_SN0_atm_ki" <- function() {
    eq_RN0 <- function(ki, Q=sqrt(240)) {
        kQ = abs(ki/Q)
        L1 = Q*ecld.ogf_star_analytic(ecld(1), kQ)
        L4 = ecld.ogf_star_analytic(ecld(4),ki)
        L1-L4
    }
    rs <- uniroot(eq_RN0, lower=-30, upper=0)
    return(rs$root) # ATM ki
}
### <---------------------------------------------------------------------->
#' @rdname ecld.quartic_SN0_atm_ki
"ecld.quartic_SN0_rho_stdev" <- function() {
    -ecld.quartic_SN0_atm_ki()/sqrt(120)
}
### <---------------------------------------------------------------------->
#' @rdname ecld.quartic_SN0_atm_ki
"ecld.quartic_SN0_skew" <- function() {
    atm_ki <- ecld.quartic_SN0_atm_ki()
    X = sqrt(abs(atm_ki))
    Z = abs(atm_ki/sqrt(240))
    sqrt(pi)*exp(Z^2-X) * (1+X-exp(X)*ecd.erfc(Z))
}
### <---------------------------------------------------------------------->
#' @rdname ecld.quartic_SN0_atm_ki
"ecld.quartic_SN0_max_RNV" <- function(sigma=0) {
    
    find_mu_p_ratio <- function(mu_p_ratio) {
        
        find_Q_SN0 <- function(Q, ki) {
            ki1 <- ki/Q - mu_p_ratio*sqrt(120)/Q # (k-mu_D1)/sigma1
            Lc1 <- ecld.ogf_star_analytic(ecld(1), ki1)
            Lc  <- ecld.ogf_star_analytic(ecld(4), ki)
            ecd.mp2f(Q*Lc1 - Lc)
        }
        
        find_skew <- function(ki, Q=sqrt(240)) {
            X = sqrt(abs(ki))
            Z = abs(ki/Q - mu_p_ratio*sqrt(120)/Q)
            sqrt(pi)*exp(Z^2-X) * (1+X-exp(X)*ecd.erfc(Z))
        }
        
        atm_ki <- uniroot(find_skew, lower=-4, upper=-0.1)$root
        Q <- uniroot(find_Q_SN0, ki=atm_ki, lower=0.1, upper=100)$root
        (Q-sqrt(240))*1000
    }
    
    if (sigma==0) {
        R <- uniroot(find_mu_p_ratio, lower=0.29, upper=0.3)$root
        return (R)
    }
    stop("ERROR!")
}
### <---------------------------------------------------------------------->
