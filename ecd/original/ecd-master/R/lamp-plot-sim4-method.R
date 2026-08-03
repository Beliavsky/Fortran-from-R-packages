#' Plot the simulation result in standard layout
#' 
#' Plot the simulation result in standard layout, with 4 or 6 charts
#' The PDF and log(PDF) histogram of Z, the lambda process.
#' The log(PDF) histogram of N, the stable count process.
#' The log(PDF) histogram of B, the binomial random walk process.
#' The 6-chart plot also includes the asymptotic kurtosis and stdev 
#' vs the bps of data points dropped in Z.
#'
#' @param object an object of lamp class
#'
#' @return an object of lamp class
#'
#' @keywords plot
#'
#' @author Stephen H-T. Lihn
#'
#' @importFrom stats sd
#' @importFrom moments kurtosis
#' 
#' @export lamp.plot_sim4
#' @export lamp.plot_sim6
#'
### <======================================================================>
"lamp.plot_sim4" <- function(object)
{
    if (length(object@Z)==0) stop("lamp object has no simulation result")
    # -----------
    layout(matrix(c(1,1,2,3,4,5),ncol=2, byrow=TRUE), 
           heights=c(1,6,6))
    
    par(mar=c(1,1,1,1))
    plot.new()
    text(0.5,0.5,lamp.plot_sim_title(object),cex=1.4,font=2)
    
    par(mar=c(5,4,4,2))
    
    lamp.plot_pdf_Z(object)
    lamp.plot_logpdf_Z(object)
    lamp.plot_logpdf_N(object)
    lamp.plot_logpdf_B(object)

}
### <---------------------------------------------------------------------->
#' @rdname lamp.plot_sim4
"lamp.plot_sim6" <- function(object)
{
    if (length(object@Z)==0) stop("lamp object has no simulation result")
    # -----------
    layout(matrix(c(1,1,2,3,4,5,6,7),ncol=2, byrow=TRUE), 
           heights=c(1,6,4,4))
    
    par(mar=c(1,1,1,1))
    plot.new()
    text(0.5,0.5,lamp.plot_sim_title(object),cex=1.4,font=2)
    
    par(mar=c(5,4,4,2))
    
    lamp.plot_pdf_Z(object)
    lamp.plot_logpdf_Z(object)
    lamp.plot_asymp_kurt(object)
    lamp.plot_asymp_sd(object)
    lamp.plot_logpdf_N(object)
    lamp.plot_logpdf_B(object)
    
}
### <---------------------------------------------------------------------->
lamp.plot_sim_title <- function(object) {
    
    lambda <- object@lambda
    beta <- object@beta
    ld0 <- ecld(lambda)
    sd_lamp <- if (is.na(object@sd)) ecld.sd(ld0) else object@sd
    k_lamp <- ecld.kurtosis(ld0)
    
    b <- sprintf("%.4f", beta)
    if (beta %in% c(0,1,-1)) b <- sprintf("%.0f", beta)
    
    ttl <- sprintf("Z: Expected L=%.2f b=%s K=%.2f sd=%.3f T.inf=%.0f wlk=%.0f", 
                      lambda, b, k_lamp, sd_lamp, object@T.inf, object@rnd.walk)
    return(ttl)
}    
### <---------------------------------------------------------------------->
"lamp.plot_pdf_Z" <- function(object)
{
    lambda <- object@lambda
    Z <- object@Z
    n <- length(Z)
    if (n==0) stop("lamp object has no simulation result in Z")
    
    sd0 <- ecld.sd(ecld(lambda))
    sd_exp <- if (is.na(object@sd)) sd0 else object@sd
    # ------------
    s <- stats::sd(Z)
    k <- moments::kurtosis(Z)
    
    lambda_z <- 2
    if (k < 100) {
        f <- function(lambda) ecld.kurtosis(ecld(lambda))-k
        lambda_z <- uniroot(f, lower=0.5, upper=6)$root
    }
    ld  <- ecld.from_sd(lambda_z, sd=s) # actual
    ld0 <- ecld.from_sd(lambda,   sd=s) # expected
    
    bins <- if (n<4000) 400 else 800
    h <- hist(Z, plot=FALSE, breaks=bins)
    # draw theoretical log-PDF from actual and expected ld
    p <- seq(-s, s, length.out=1000)
    y1 <- ecld.pdf(ld, p)
    y2 <- ecld.pdf(ld0, p)
    
    ymax <- max(c(h$density,y1,y2), na.rm=TRUE)
    plot(h$mids, h$density, col="blue", cex=0.3,
         xlim = c(-s, s), 
         ylim = c(0, ymax),
         xlab="Z",
         ylab="PDF(Z)",
         main=sprintf("Z: rlz L=%.2f K=%.2f sd=%.5f bins=%d", 
                      lambda_z, k, s, bins))
    
    lines(p, y1,  col="red")
    lines(p, y2, col="purple")
    
    abline(v=0, lty=2)
    abline(v=s*0.5, lty=2)
    abline(v=s, lty=2)
    
    text(s*0.5, ymax*0.8, pos=2, label=sprintf("sd/2"))
    text(s, ymax*0.8, pos=2, label=sprintf("sd"))
    
}
### <---------------------------------------------------------------------->
"lamp.plot_logpdf_Z" <- function(object)
{
    lambda <- object@lambda
    beta <- object@beta
    Z <- object@Z
    n <- length(Z)
    if (n==0) stop("lamp object has no simulation result in Z")
    
    sd0 <- ecld.sd(ecld(lambda))
    sd_exp <- if (is.na(object@sd)) sd0 else object@sd
    # ------------
    s <- stats::sd(Z)
    k <- moments::kurtosis(Z)
    
    lambda_z <- 2
    if (k < 100) {
        f <- function(lambda) ecld.kurtosis(ecld(lambda))-k
        lambda_z <- uniroot(f, lower=0.5, upper=6)$root
    }
    ld  <- ecld.from_sd(lambda_z, sd=s) # actual
    ld0 <- ecld.from_sd(lambda,   sd=s) # expected

    bins <- if (n<4000) 200 else 400
    h <- hist(Z, plot=FALSE, breaks=bins)
    
    # LM to figure out lambda
    x <- h$mids
    y <- log(h$density)
    
    lx1 <- log(abs(x))
    ly1 <- log(-y)
    I <- which(is.finite(lx1) & is.finite(ly1))
    lx2 <- lx1[I]
    ly2 <- ly1[I]
    
    slm <- lm(ly2 ~ lx2)
    clm <- slm$coefficients
    lambda_imp <- 2/clm[2]
    
    #plot(lx2, ly2, cex=0.3)
    #p <- seq(min(lx2), max(lx2), length.out=100)
    # lines(p, clm[1]+p*clm[2], col="red")
    
    ymax <- max(y, na.rm=TRUE)
    plot(x, y, col="blue", type="l", cex=0.3,
         xlab="Z",
         ylab="log(PDF(Z))",
         main=sprintf("Z: rlz L=%.2f K=%.2f sd=%.5f bins=%d", 
                      lambda_z, k, s, bins))
    
    # draw theoretical log-PDF from actual and expected ld
    p <- seq(min(h$mids), max(h$mids), length.out=1000)
    lines(p, log(ecld.pdf(ld, p)),  col="red")
    lines(p, log(ecld.pdf(ld0, p)), col="purple")
    
    abline(v=0, lty=2)
    abline(v=s*10, lty=2)
    text(s*10, ymax-2, pos=2, label=sprintf("10 sd"))
    
    if (abs(beta) > 0 & beta < 1) {
        text(0, ymax-0.5, pos=4, 
             label=sprintf("    skewness %.3f", moments::skewness(Z)))
    }
    # this doesn't appear to be helpful
    #text(0, max(y, na.rm=TRUE)*1.25, pos=4, cex=1,
    #     label=sprintf("    slope implied lambda %.2f K=%.2f", 
    #                   lambda_imp, ecld.kurtosis(ecld(lambda_imp))))
}
### <---------------------------------------------------------------------->
"lamp.plot_logpdf_N" <- function(object)
{
    N <- object@N
    lambda <- object@lambda
    if (length(N)==0) stop("lamp object has no simulation result in N")
    
    # ------------
    h_n <- hist(N, breaks=200, plot=FALSE)
    # this may have issues for scaling? need a better logic
    dscum <- cumsum(h_n$density)*diff(h_n$mids)[1]
    pct <- 0.97 # ignore the tail for the lm fit
    I <- which(is.finite(log(h_n$density)) & h_n$mids > 0.5 & dscum < pct)
    c <- s <- 1
    n_max_lm <- max(h_n$mids[I], na.rm=TRUE)
    # this only works for L=4 in real world !?
    if (lambda > 2) {
        h_n_lm <- lm(log(h_n$density)[I] ~ h_n$mids[I], na.action = na.omit) 
        c <- h_n_lm$coefficients[1]
        s <- h_n_lm$coefficients[2]
    }
    plot(h_n$mids, log(h_n$density), col="blue", type="l",
         xlab="N",
         ylab="log(PDF(N))",
         main=sprintf("N: sim=%d slope=%.3f mean=%.2f", length(N), s, mean(N)))
    lines(h_n$mids, c+s*h_n$mids, col="red")
    abline(v=n_max_lm, lty=2)
    abline(v=0, lty=2)
    abline(v=object@N.lower, lty=2)
    text(n_max_lm, max(log(h_n$density)[I], na.rm=TRUE)-1, pos=4, 
         label=sprintf("%.2f pct", pct*100))
}
### <---------------------------------------------------------------------->
"lamp.plot_logpdf_B" <- function(object)
{
    B <- object@B
    lambda <- object@lambda
    beta <- object@beta
    if (length(B)==0) stop("lamp object has no simulation result in B")
    # ------------
    h_b <- hist(B, breaks=200, plot=FALSE)
    s <- stats::sd(B)
    k <- moments::kurtosis(B)
    lambda_b <- 1
    if (k > 3 & k < 10) {
        f <- function(lambda) ecld.kurtosis(ecld(lambda))-k 
        # no analytic formula for skewed kurtosis! beta=beta_b
        lambda_b <- uniroot(f, lower=0.5, upper=11)$root
    }
    m <- max(log(h_b$density), na.rm=TRUE)
    ld <- ecld.from_sd(lambda=lambda_b, sd=s)
    ld1 <- ld
    if (object@beta==0.005) {
        beta_b <- 0.21
        mu_b <- 0.1
        ld1@beta <- beta_b # again unfortunately, skew is not defined for lambda < 2
        ld1@mu <- mu_b
        ld1@is.sged <- TRUE
    }
    ymax <- max(log(h_b$density), na.rm=TRUE)
    plot(h_b$mids, log(h_b$density), col="blue", type="l",
         xlab="B",
         ylab="log(PDF(B))",
         main=sprintf("B: L=%.3f sd=%.5f K=%.2f wlk=%.0f", 
                      lambda_b, s, k, object@rnd.walk))
    lines(h_b$mids, log(ecld.pdf(ld1, h_b$mids)), col="red")
    abline(v=0, lty=2)
    
    if (abs(beta) > 0 & beta < 1) {
        text(0, ymax-0.5, pos=4, 
             label=sprintf("    skewness %.3f", moments::skewness(B)))
    }
}
### <---------------------------------------------------------------------->
"lamp.plot_asymp_kurt" <- function(object)
{
    Z <- object@Z
    lambda <- object@lambda
    n <- length(Z)
    if (n==0) stop("lamp object has no simulation result in Z")
    
    calc_kurtosis <- function(drop, Z) {
        Z1 <- lamp.drop_outliers(Z, drop)
        moments::kurtosis(Z1)
    }
    # ------------
    ld0 <- ecld(object@lambda)
    k_exp <- ecld.kurtosis(ld0)

    samples <- if (n >= 100) 0:100 else 0:n
    KK <- ecd.mcsapply(ld0, samples, function(d) calc_kurtosis(d,Z))
    bps <- (1:length(KK))/n*10000
    ymin <- min(KK)
    ymax <- max(KK)
    if (ymax > k_exp*0.8) ymax <- max(c(ymax, k_exp))
    plot(bps, KK, cex=0.3, col="blue",
         ylim=c(ymin,ymax),
         xlab="bps of drops",
         ylab="kurtosis(Z)",
         main=sprintf("Z: kurt exp'ed %.2f rlz %.2f", 
                      k_exp, max(KK)))
    abline(h=k_exp, lty=2)
    
}
### <---------------------------------------------------------------------->
"lamp.plot_asymp_sd" <- function(object)
{
    Z <- object@Z
    lambda <- object@lambda
    n <- length(Z)
    if (n==0) stop("lamp object has no simulation result in Z")
    
    calc_sd <- function(drop, Z) {
        Z1 <- lamp.drop_outliers(Z, drop)
        stats::sd(Z1)
    }
    # ------------
    ld0 <- ecld(object@lambda)
    sd_lamp <- if (is.na(object@sd)) ecld.sd(ld0) else object@sd
    
    samples <- if (n >= 200) 0:200 else 0:n
    SS <- ecd.mcsapply(ld0, samples, function(d) calc_sd(d,Z))
    bps <- (1:length(SS))/n*10000
    ymin <- min(SS)
    ymax <- max(SS)
    if (ymax > sd_lamp*0.8) ymax <- max(c(ymax, sd_lamp))
    plot(bps, SS, cex=0.3, col="blue",
         ylim=c(ymin,ymax),
         xlab="bps of drops",
         ylab="sd(Z)",
         main=sprintf("Z: sd exp'ed %.3f rlz %.3f", 
                      sd_lamp, max(SS)))
    abline(h=sd_lamp, lty=2)
    
}
### <---------------------------------------------------------------------->
lamp.drop_outliers <- function(x, drop=1)
{
    if (drop==0) return(x)
    if (drop<0) stop("Drop must be a positive integer")
    
    `%notin%` <- function (x, table) match(x, table, nomatch = 0L) == 0L
    z <- abs(x)
    z1 <- rev(z[order(z)])
    x[abs(x) %notin% utils::head(z1, drop)]
}
### <---------------------------------------------------------------------->
