#' Calculate ATM attributes from key quartic parameters
#'
#' This utility takes a data frame of key quartic parameters, and generates several key ATM attributes.
#' Input fields are: ttm - time to expiration, sigma - term structure of sigma,
#' epsilon_ratio - term structure of epsilon/sigma, mu_plus_ratio - term structure of (mu_p-mu_D)/stdev.
#' The output fields are: atm_ki, atm_kew, atm_vol, rho, and rho_ratio - rho/stdev.
#'
#' @param df data.frame
#' @param dt character, one of three sample dates used in the quartic model paper (YYYY-MM-DD)
#' @param ttm numeric, list of time to expiration (T=1 for one year)
#' @param target_file character, file location to cache the attribute data (to avoid lengthy repetitions)
#' @param skew_adjusted logical, if true, use skew adjusted T=0 intercep,
#'                      else use the tercep from linear fit. Default is \code{TRUE}.
#'
#' @return data.frame
#'
#' @keywords quartic
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecld.quartic_Qp_atm_attr
#' @export ecld.quartic_model_sample
#' @export ecld.quartic_model_sample_attr
#'
#' @importFrom utils read.table
#' @importFrom utils write.table
#'
#' @examples
#' ttm <- seq(sqrt(90), sqrt(365), length.out=3)^2 / 365
#' epsr = 0.014 + 0*ttm
#' mupr <- -(ecld.quartic_SN0_max_RNV() - 0.2*sqrt(ttm))
#' \dontrun{
#'     df <- data.frame(ttm=ttm, sigma=0.2*sqrt(ttm/120), mu_plus_ratio=mupr, epsilon_ratio=epsr)
#'     ecld.quartic_Qp_atm_attr(df)
#' }
### <======================================================================>
"ecld.quartic_Qp_atm_attr" <- function(df) {
    
    if (! is(df, "data.frame")) stop("Error: expect df to be a data.frame")
    if (nrow(df) != 1) {
        dfl <- split(df, 1:nrow(df))
        df2 <- do.call("rbind", parallel::mclapply(dfl, ecld.quartic_Qp_atm_attr))
        return (df2)
    }
    
    ttm <- df$ttm
    sigma <- df$sigma
    epsilon_ratio <- df$epsilon_ratio
    mu_plus_ratio <- df$mu_plus_ratio
    
    df$atm_ki <- NaN
    df$atm_skew <- NaN
    df$atm_vol <- NaN
    df$rho_ratio <- NaN
    df$rho <- NaN
    
    ld0 <- tryCatch(
        ecld.quartic(sigma, epsilon_ratio*sigma, rho=0, mu_plus_ratio),
        error = function(cond) { return(NULL) }
    )
    if (is.null(ld0)) return(df)
    
    df$atm_ki <- tryCatch(
        ecld.quartic_Qp_atm_ki(ld0),
        error = function(cond) { return(NaN) }
    )
    if (is.nan(df$atm_ki)) return(df)

    df$atm_skew <- tryCatch(
        ecld.quartic_Qp_skew(ld0, df$atm_ki)/sqrt(2*ttm),
        error = function(cond) { return(NaN) }
    )
    if (is.nan(df$atm_skew)) return(df)
    
    ld <- ld0
    rho <- tryCatch(
        ecld.quartic_Qp_rho(ld0, df$atm_ki),
        error = function(cond) { return(NaN) }
    )
    if (is.nan(rho)) return(df)

    ld@rho <- rho
    df$rho_ratio <- ecd.mp2f(ld@rho/ecld.sd(ld))
    df$rho <- ld@rho

    df$atm_vol <- tryCatch(
        ecd.mp2f(ecld.op_VL_quartic(ld, 0, otype="p", ttm=ttm)),
        error = function(cond) { return(NaN) }
    )

    df    
}
### <---------------------------------------------------------------------->
#' @rdname ecld.quartic_Qp_atm_attr
"ecld.quartic_model_sample" <- function(dt, ttm, skew_adjusted=TRUE) {
    
    IV_atm <- NULL
    S_d_mu <- NULL
    epsilon_ratio <- NULL
    if (dt == "2015-07-20") {
        IV_atm <- 0.08*(1+sqrt(ttm))
        epsilon_ratio = 0.0136
        mu_plus_ratio <- 0.29 -0.479*ttm^0.5
        if (skew_adjusted) {
            mu_plus_ratio <- 0.230 -0.53*ttm^0.5
            # mu_plus_ratio <- 0.220 -0.31*ttm # 0.8*ttm^0.75
            # mu_plus_ratio <- ifelse(mu_plus_ratio < 0, 8*mu_plus_ratio, mu_plus_ratio)
        }
    }
    else if (dt == "2015-08-24") {
        IV_atm <- 0.2*ttm^-0.27
        epsilon_ratio = 0.0174
        mu_plus_ratio <- 0.295 -0.185*ttm^0.5
        if (skew_adjusted) {
            mu_plus_ratio <- 0.315 -0.21*ttm^0.5
            # mu_plus_ratio <- 0.30 -0.24*ttm
        }
    }
    else if (dt == "2015-10-06") {
        IV_atm <- 0.18
        epsilon_ratio = 0.0094
        mu_plus_ratio <- 0.293 -0.337*ttm^0.5
        if (skew_adjusted) {
            mu_plus_ratio <- 0.325 -0.44*ttm^0.5
            # mu_plus_ratio <- 0.293 -0.45*ttm
            # mu_plus_ratio <- ifelse(mu_plus_ratio < 0, 8*mu_plus_ratio, mu_plus_ratio)
        }
    }
    else {
        stop(paste("Error: date not in scope:", dt))
    }
    
    df <- data.frame(ttm=ttm,
                     sigma=IV_atm*sqrt(ttm/120),
                     mu_plus_ratio=-mu_plus_ratio, # negative sign for puts
                     epsilon_ratio=epsilon_ratio)
    df
}
### <---------------------------------------------------------------------->
#' @rdname ecld.quartic_Qp_atm_attr
"ecld.quartic_model_sample_attr" <- function(dt, ttm, target_file, skew_adjusted=TRUE) {
    
    if (is.null(ttm)) ttm <- seq(sqrt(0.5), sqrt(365), length.out=40)^2 / 365

    df <- ecld.quartic_model_sample(dt, ttm, skew_adjusted)
    md5 <- digest::digest(df, "md5")
    
    if (file.exists(target_file)) {
        df2 <- utils::read.table(file=target_file)
        md5_2 <- unique(df2$md5)
        if (length(md5_2)==1 & md5==md5_2) return(df2)
    }
    df$md5 <- md5
    df <- ecld.quartic_Qp_atm_attr(df)
    utils::write.table(df, target_file)
    df
}
### <---------------------------------------------------------------------->
