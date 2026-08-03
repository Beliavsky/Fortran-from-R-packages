#' Produce 3x3 plot of VIX volatility smiles for a date
#'
#' This utlity produces 3x3 plot of volatility smiles for a date. It is used for the VIX option paper.
#'
#' @param option_data dataframe, read from \code{ecop.read_csv_by_symbol}
#' @param result dataframe, the VIX optimx result
#' @param result_avg dataframe, the VIX optimx result using average lambda for all expirations
#' @param date_str character in the form of YYYY-MM-DD
#'
#' @return The 3x3 plot
#'
#' @keywords term structure
#'
#' @export ecop.vix_plot_3x3
#'
### <======================================================================>
"ecop.vix_plot_3x3" <- function(date_str, option_data, result, result_avg)
{
    # for R CMD check - no visible binding for global variable
    CLS <- NULL
    Date <- NULL
    PC <- NULL
    datadt <- NULL
    
    idx_range <- 1:8
    df <- subset(option_data, Date == as.Date(date_str) & CLS == "VIX")
    days <- sort(unique(df[["days"]]))
    vols <- sapply(days, function(d) sum(subset(df, days==d)[["VOL"]])/1000.0)
    
    days <- days[vols>0] # remove days that have no volume!
    vols <- vols[vols>0]

    get_S <- function(d) {
        df_c <- subset(df, days==d & PC=="C")
        head(unique(df_c$UNDL_PRC),1)
    }
    
    dtn <- as.numeric(format(as.Date(date_str), "%Y%m%d"))
    
    for (idx in idx_range) {
        
        i <- idx
        rs <- subset(result, datadt == dtn & idx == i)
        rs_avg <- subset(result_avg, datadt == dtn & idx == i)
        ttm <- rs$ttm
        
        d <- days[idx]

        min_k = -0.8
        max_k = 0.8
        
        df_c <- subset(df, days==d & PC=="C")
        df_p <- subset(df, days==d & PC=="P")
        S <- S_raw <- head(unique(df_c$UNDL_PRC),1)

        K0_c <- df_c$STRK_PRC
        K0_p <- df_p$STRK_PRC
        IV_c <- df_c$IVOL
        IV_p <- df_p$IVOL

        k_p <- log(K0_p/S)
        valid_p <- which(k_p >= min_k & k_p <= max_k & IV_p > 0)
        IV_p <- IV_p[valid_p]
        k_p <- k_p[valid_p]
        
        k_c <- log(K0_c/S)
        valid_c <- which(k_c >= min_k & k_c <= max_k & IV_c > 0)
        IV_c <- IV_c[valid_c]
        k_c <- k_c[valid_c]
        
        # ---------------------------------------------------------------
        otype <- "c"
        beta <- 0
        lambda <- rs$lambda
        mu_plus_ratio <- rs$mu_plus_ratio_r/100
        epsilon_ratio <- rs$epsilon_ratio_r/100
        sd <- rs$sd
        rho <- rs$rho
        
        # ----------------------------------------------------------------
        ld1 <- ecop.get_ld_triple(lambda, sd, beta, mu_plus_ratio, epsilon_ratio, atm_imp_k=0)$ld1
        ld1@rho <- rho
        k_all = sort(unique(c(k_p, k_c)))
        L0 <- ecd.mp2f(ecld.ogf(ld1, k_all-rho, otype=otype, RN=FALSE))
        L <- ecd.mp2f(L0 + ld1@epsilon)
        
        ivol_ecop  <- 1/sqrt(2*ttm) * ecd.mp2f(ecld.op_V(L,  k_all-rho, otype=otype))
        ivol0_ecop <- 1/sqrt(2*ttm) * ecd.mp2f(ecld.op_V(L0, k_all-rho, otype=otype))

        # experimental, using avg_lambda for all curves
        avg_lambda <- rs_avg$lambda

        ld2 <- ecop.get_ld_triple(avg_lambda, rs_avg$sd, beta,
                mu_plus_ratio = rs_avg$mu_plus_ratio_r/100,
                epsilon_ratio = rs_avg$epsilon_ratio_r/100,
                atm_imp_k=0)$ld1
        ld2@rho <- rs_avg$rho
        L0_2 <- ecd.mp2f(ecld.ogf(ld2, k_all-ld2@rho, otype=otype, RN=FALSE))
        L2 <- ecd.mp2f(L0_2 + ld2@epsilon)

        ivol2_ecop  <- 1/sqrt(2*ttm) * ecd.mp2f(ecld.op_V(L2,  k_all-ld2@rho, otype=otype))

        # ---------------------------------------------------------------
        IV_max <- max(c(IV_p, IV_c, ivol_ecop), na.rm=TRUE)
        IV_min <- min(c(IV_p, IV_c, ivol_ecop), na.rm=TRUE)

        IV_ptr <- mean(ivol_ecop[k_all < 0], na.rm=TRUE)

        min_k <- min(k_all) # for plot
        # ---------------------------------------------------------------
        
        plot(k_p, IV_p, type="p", pch=1, col="blue", lwd=2, cex=0.7,
             xlab="$k$ (log-strike)", ylab="$\\sigma_{BS}(k)$",
             xlim=c(min_k, max_k),
             ylim=c(IV_min*0.8, IV_max*1.1),
             main=sprintf("%dd, T=%.3f", d, ttm))
        
        
        points(k_c, IV_c, col="green", cex=0.7)
        lines(k_all, ivol_ecop, col="red", lwd=3)
        lines(k_all, ivol2_ecop, col="orange", lwd=2, lty=2) # experimental
        
        abline(v=0, col="black", lty=2)
        segments(rho, 0, rho, IV_ptr)
        segments(rho+sd, 0, rho+sd, IV_ptr)
        segments(0, sd/sqrt(ttm), 0.2, sd/sqrt(ttm))

        text(x=-0.5, y=IV_max*0.9, labels=
            sprintf("$\\lambda$=%.2f/%.2f", lambda, avg_lambda), cex=0.8, pos=4)


    }
    
    plot(days, sapply(days, get_S), col="blue", lwd=2, cex=0.7,
         xlab="days to expire", ylab="px",
         main=sprintf("VIX futures %s", date_str))
    
}
### <---------------------------------------------------------------------->
