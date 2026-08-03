#' Produce 3x3 plot of volatility smiles for a date
#'
#' This utlity produces 3x3 plot of volatility smiles for a date. It is used for the term structure paper.
#'
#' @param term_data term structure data for one date, produced from \code{ecop.term_master_calculator}
#' @param date_str character in the form of YYYY-MM-DD
#' @param trim_points integer, specifying number of data points to present in the plots
#' @param days list of days to expiration from market data
#' @param target_days list of ceiling days for the plot
#' @param realized_days list of days realized for the plot
#' @param add.first.day logic, whether to add the first expiration date to \code{target_days}.
#'                      Default is \code{TRUE}.
#' @param show.put.bid logic, show bid smile for put option. Default is \code{FALSE}.
#' 
#' @return The 3x3 plot
#'
#' @keywords term structure
#'
#' @export ecop.term_plot_3x3
#' @export ecop.term_target_days_default 
#' @export ecop.term_realized_days 
#' @export ecop.term_idx_range
#'
### <======================================================================>
"ecop.term_plot_3x3" <- function(term_data, date_str, trim_points=151, 
                                 target_days=NULL, add.first.day=TRUE, show.put.bid=FALSE)
{
    
    config <- NULL
    if ("quartic.config" %in% names(term_data)) config <- term_data$quartic.config
    get_config_by_idx <- function(i) { subset(config, idx==i) }
    
    days <- term_data$days

    if (is.null(target_days)) target_days <- ecop.term_target_days_default
    realized_days <- ecop.term_realized_days(target_days, days)
    idx_range <- ecop.term_idx_range(realized_days, days)
    if (add.first.day) idx_range <- c(1,idx_range)
    
    for (idx in idx_range) { # up to 27
        
        put <- term_data$dt[[idx]]$put
        call <- term_data$dt[[idx]]$call
        d <- days[idx]
        
        min_k = -0.6
        max_k = 0.3
        if (d > 120) min_k <- -0.8
        
        k_p <- put$k
        valid_p <- which(k_p >= min_k & k_p <= max_k)
        IV_p <- put$IV[valid_p]
        k_p <- k_p[valid_p]
 
        k2_p <- put$k2
        valid2_p <- which(k2_p >= min_k & k2_p <= max_k & put$IV2 > 0)
        IV2_p <- put$IV2[valid2_p]
        k2_p <- k2_p[valid2_p]
        
        k_c <- call$k
        valid_c <- which(k_c >= min_k & k_c <= max_k)
        IV_c <- call$IV[valid_c]
        k_c <- k_c[valid_c]

        k2_c <- call$k2
        valid2_c <- which(k2_c >= min_k & k2_c <= max_k & call$IV2 > 0)
        IV2_c <- call$IV2[valid2_c]
        k2_c <- k2_c[valid2_c]
        
        # ---------------------------------------------------------------
        # trim data to make plots look cleaner
        N = trim_points
       
        trim_data <- function(N, k, IV, IV2) {
            min_k <- min(k)
            max_k <- max(k)
            
            # k1 <- k[abs(k) <= 0.1] # includes all points with 10%
            k2 <- seq(min_k, max_k, length.out=N)
            k3 <- c()
            for (i in 1:N) {
                k3[i] <- min(k[k >= k2[i]])
            }
            # k3 = sort(unique(c(k1,k3)))
            k3 = sort(unique(k3))
            
            IV <- IV[k %in% k3]
            IV2 <- IV2[k %in% k3]
            k <- k3
            return(list(k=k, IV=IV, IV2=IV2))
        }
        
        trim_p <- trim_data(N, k_p, IV_p, IV2_p)
        k_p <- trim_p$k
        IV_p <- trim_p$IV
        IV2_p <- trim_p$IV2
        
        trim_c <- trim_data(N, k_c, IV_c, IV2_c)
        k_c <- trim_c$k
        IV_c <- trim_c$IV
        IV2_c <- trim_c$IV2
        
        # ---------------------------------------------------------------
        
        IV_max <- max(c(IV_p), na.rm=TRUE)
        IV_min <- min(c(IV_p, IV_c, IV2_p[IV2_p > 0.02]), na.rm=TRUE)
        # ---------------------------------------------------------------
 
        pm = get_config_by_idx(idx)
        sigma = pm[,"sigma_r"] * 1e-3 * ecd.mp1
        mu_plus_ratio = pm[,"mu_plus_ratio_r"] * (-1) # puts
        epsilon = pm[,"epsilon_r"] * 1e-5
        rho = pm[,"r_m_r"] 
        ttm = pm[,"ttm"] 
        
        # ---------------------------------------------------------------
        
        plot(k_p, IV_p, type="p", pch=1, col="blue", lwd=2, cex=0.7,
             xlab="$k$ (log-strike)", ylab="$\\sigma_{BS}(k)$",
             xlim=c(min_k, max_k), 
             ylim=c(IV_min, IV_max),
             main=sprintf("%dd, T=%.3f", d, ttm))
        
        if (show.put.bid) points(k2_p, IV2_p, col="orange", pch=3)
        
        lines(k_c, IV_c, col="green", lwd=2)
        points(k_p, IV2_p, col="orange", cex=0.3)
        
        abline(v=0, col="black", lty=2)

        # ----------------------------------------------------------------

        ld1 <- ecld.quartic(sigma=sigma, rho=rho, epsilon=epsilon, mu_plus_ratio=mu_plus_ratio)
        k_all = sort(unique(c(k_p, k_c)))
        ivol_ecop <- ecd.mp2f(ecld.op_VL_quartic(ld1, k_all, otype="p", ttm=ttm))
        
        lines(k_all, ivol_ecop, col="red", lwd=3)
        
        y_sd <- ecd.mp2f(ecld.sd(ld1))/sqrt(ttm)
        
        segments(-0.1, y_sd, 0, y_sd)
        text(x=-0.1, y=y_sd, labels=sprintf("stdev $\\sqrt{T}$ = %.3f", y_sd), cex=0.8, pos=2)

    }
    
}
### <---------------------------------------------------------------------->
#' @rdname ecop.term_plot_3x3
ecop.term_target_days_default <- c(14,35, 65,100,130, 160,270,365)
### <---------------------------------------------------------------------->
#' @rdname ecop.term_plot_3x3
ecop.term_realized_days <- function(target_days, days) {
    sapply(target_days, function(d) max(days[days <=d]))
}
### <---------------------------------------------------------------------->
#' @rdname ecop.term_plot_3x3
ecop.term_idx_range <- function(realized_days, days) {
    which(days %in% realized_days)
}
### <---------------------------------------------------------------------->
