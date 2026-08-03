#' Fit observations to QSLD via MLE
#' 
#' This utility will fit the observations to qsld with MLE using \code{optimx}. 
#' There are three features:
#' First, it has the ability to provide initial estimate to save the user 
#' from the headache of guessing. 
#' Second, the user has the flexibility of fixing convolution, \code{beta.a}, and/or \code{nu0/theta} ratio.
#' And the user can ask the utility to estimate mu as well.
#' Third, the MLE optimization comes with two flavors: 
#' (a) the log-likelihood calculated from the PDF of all observations; or
#' (b) the log-likelihood calculated from the histogram of the observations.
#' The later is much faster than the former. One can use the later to obtain a good estimate,
#' then feed it into the former, if necessary. 
#' This utility also comes with an LMS regression method called "pdf.lms". This option can be used 
#' as a preprocessor to the MLE methods. It regresses the theoretical PDF against the empirical PDF
#' obtained from the histogram, and minimizes the LMS of PDF difference within 2-stdev 
#' and log(PDF) difference within 4-stdev.
#'
#' @param x numeric, the observation of log-returns.
#' @param breaks numeric, the breaks for the histogram of observations. 
#'               For \code{lme}, this parameter is only for display purpose.
#' @param init.qsld an object of sld class as an initial qsld guess. The user can request the utility to 
#'                  estimate the initial parameters by setting nu0 to \code{NaN}. However, t must be provided.
#' @param method character, optimization algorithm to use, it could be either
#'               \code{fast.mle}, \code{mle}, or \code{pdf.lms}. Default is \code{fast.mle}.
#' @param fix.convo     numeric, fix convolution to a specific number, default is \code{NaN}.
#' @param fix.beta.a    numeric, fix annualized beta to a specific number, default is \code{NaN}.
#' @param fix.nu0.ratio numeric, fix \eqn{\nu_0/\theta} ratio to a specific number, default is \code{NaN}.
#' @param derive.mu     logical, if specified, to derive \code{mu} automatically, default is \code{TRUE}.
#' @param plot.interval numeric, interval of iterations to plot the fit, default is 20.
#'                      If set to zero, the plot is disabled.
#' @param verbose.interval numeric, interval of iterations to print verbose message, default is 20.
#'                      If set to zero, the verbose message is disabled.
#' @param itnmax    numeric, specify maximum iterations for optimx, default is 500. 
#'
#' @return a list of two components: \code{qsld} as an object of sld class representing the QSLD fit;
#'         \code{optimx.out} storing the raw output from \code{optimx}.
#'
#' @keywords fit
#'
#' @author Stephen H-T. Lihn
#'
#' @importFrom optimx optimx
#' @importFrom moments kurtosis
#' 
#' @export qsld.fit
#'
#' @examples
#' \dontrun{
#'   x <- ecd.data.arr("spx", lag=1, drop=1)$x
#'   breaks <- 200
#'   t <- 1/250
#'   d <- qsld(t=t)
#'   d@nu0 <- NaN # request utility to estimate
#    output <- qsld.fit(x, breaks=breaks, init.qsld=d, derive.mu=TRUE)
#' }
### <======================================================================>
"qsld.fit" <- function(x, breaks, init.qsld, method = "fast.mle",
                       fix.convo = NaN, fix.beta.a = NaN, fix.nu0.ratio = NaN, derive.mu = TRUE,
                       plot.interval = 20, verbose.interval = 20, itnmax = 500) {

    d <- init.qsld
    t <- d@t
    h <- hist(x, breaks=breaks, plot=FALSE)
    dx0 <- diff(h$mids)[1]
    mean_x <- mean(x)
    sd_x <- sd(x)
    length_x <- length(x)
    kurt_x <- moments::kurtosis(x)
    
    # ---------------------------------------------------------------------
    # if nu0 is NaN, then we estimate the guess
    if (is.na(d@nu0)) {
        pdf0 <- max(hist(x,breaks=breaks, plot=FALSE)$density, na.rm=TRUE) * sd(x)
        nr <- if (is.na(fix.nu0.ratio)) 6 else fix.nu0.ratio
        theta <- sqrt(var(x)/t/((nr+6)^2+24)) # assume nu0 = 6 theta
        
        # estimate convo
        convo = 1
        repeat {
            pdf1 <- qsl_std_pdf0_analytic(t=t, nu0=theta*nr, theta=theta, convo=convo, beta.a=0) 
            if (verbose.interval > 0) print(sprintf("estimating theta %.2f convo %.2f pdf0 %.3f vs %.3f", theta*100, convo, pdf0, pdf1))
            if (pdf1 <= pdf0 | convo > 10) break 
            convo <- convo+0.5
        }
        d <- qsld(t=t, mu=mean(x), nu0=theta*nr, theta=theta, convo=convo, beta.a=0)
    }
    # ---------------------------------------------------------------------
    # utility to handle parameters in-and-out of optimx
    stk <- list( # a stack structure
        pm = c() # pm is input
        ,pos = rep(NaN, 5) # position of parameters, in the order of mu, nu0, theta, convo, beta.a
    )
    stk <- .stk.from_qsld(stk, d, derive.mu, fix.convo, fix.beta.a, fix.nu0.ratio)
    if (verbose.interval > 0) {    
        print("input to optimx:")
        print(stk)
        print("-------------------")
    }
    
    # --- fit ---
    iter <- 0
    fn <- function(v0) {
        iter <<- iter + 1
        #print(v0)
        v <- .stk.organize(v0, stk$pos, t, derive.mu, fix.convo, fix.beta.a, fix.nu0.ratio)
        #print(v)
        
        if (abs(v[1]) >= 1) return(1000) # mu
        if (v[2] <= 0.00001) return(2010) # nu0 
        if (v[3] <= 0.00001) return(2020) # theta 
        if (v[4] <= 0.01) return(3020) # convo
        if (abs(v[5]) >= 10) return(3010) # beta.a
        
        if (derive.mu) {
            d <- qsld(t=t, mu=0, nu0=v[2], theta=v[3], convo=v[4], beta.a=v[5])
            mean_d <- mean_x - sld.mean(d)
            v[1] <- mean_d # add mean
        }
        
        if (verbose.interval > 0 & iter %% verbose.interval == 0) {
            print(sprintf("%d mu.a %.4f nu0 %.2f theta %.3f convo %.3f beta.a %.2f", 
                      iter, v[1]/t, v[2]*100, v[3]*100, v[4], v[5]))
        }
        if (abs(v[1]) >= 1) return(1000) # annualized mu
        
        # -------------------------------------------------------------------------
        pdf_sd = 2 
        logpdf_sd = 4

        px0 <- h$mids
        p <- dqsl(px0-v[1], t=t, nu0=v[2], theta=v[3], convo=v[4], beta.a=v[5])
        
        delta <- NULL
        msg <- NULL
        if (method == "fast.mle") {
            llk <- log(p)*h$density*dx0
            llk <- llk[is.finite(llk)]
            llk <- sum(llk, na.rm=TRUE)
            delta <- -llk
            msg <- sprintf("%d %s delta: -llk= %.4f", iter, method, delta)
        }
        if (method == "mle") {
            pp <- dqsl(x-v[1], t=t, nu0=v[2], theta=v[3], convo=v[4], beta.a=v[5])
            llk <- sum(log(pp))/length_x
            delta <- -llk
            msg <- sprintf("%d %s delta: -llk= %.4f (length %d)", iter, method, delta, length_x)
        }
        # pdf.lms is a hidden method to obtain estimate based on PDF and log(PDF)
        if (method == "pdf.lms") {
            pdf_amp = 10
            
            I2 <- which(abs(h$mids-mean_x) <= pdf_sd*sd_x)
            I4 <- which(abs(h$mids-mean_x) <= logpdf_sd*sd_x)
            
            delta2 <- ((h$density-p)*sd_x*pdf_amp)[I2]
            delta2 <- delta2[is.finite(delta2)]
            delta2 <- sqrt(sum(delta2^2, na.rm=TRUE))
            
            delta4 <- (log(h$density*sd_x)-log(p*sd_x))[I4]
            delta4 <- delta4[is.finite(delta4)]
            delta4 <- sqrt(sum(delta4^2, na.rm=TRUE))
            
            delta <- delta2 + delta4 
            msg <- sprintf("%d %s delta: pdf= %.4f logpdf= %.4f total= %.4f", iter, method, delta2, delta4, delta)
        }
        if (verbose.interval > 0 & iter %% verbose.interval == 0) print(msg)
        
        if (is.null(delta)) {
            print(paste("Delta is null, unexpected error for method:", method))
            return(delta)
        }

        # plot
        if (plot.interval > 0 & iter %% plot.interval == 0) {

            px1 <- seq(-2, 2, length.out=200)*sd_x
            p1 <- dqsl(px1-v[1], t=t, nu0=v[2], theta=v[3], convo=v[4], beta.a=v[5])
            ymax <- max(c(h$density, p1), rm.na=TRUE)
            
            par(mfrow=c(1,2));
            
            plot(h$mids/sd_x, h$density*sd_x, cex=0.3, 
                 xlim=c(-1,1)*4,
                 ylim=c(0, ymax*sd_x),
                 xlab="log-return x sd", ylab="PDF/sd",
                 main=sprintf("iter %d delta %.6f", iter, delta))
            
            lines(px1/sd_x, p1*sd_x, col="orange", lty=1)
            lines(px0/sd_x, p*sd_x, col="red")
            abline(v=pdf_sd, col="orange", lty=2)
            abline(v=-pdf_sd, col="orange", lty=2)
            
            plot(h$mids/sd_x, log(h$density*sd_x), 
                 xlab="log-return x sd", ylab="log(PDF/sd)",
                 cex=0.3, main=sprintf("sd %.5f kurtosis %.2f", sd_x, kurt_x))
            lines(px1/sd_x, log(p1*sd_x), col="orange", lty=1)
            lines(px0/sd_x, log(p*sd_x), col="red")
            abline(v=logpdf_sd, col="orange", lty=2)
            abline(v=-logpdf_sd, col="orange", lty=2)
        }
        # output
        return(delta)
    }
    # ---------------------------------------------------------------------
    # print(stk$pm)
    pout <- optimx::optimx(stk$pm, fn, method = c("BFGS"), itnmax = itnmax)
    
    L <- length(stk$pm)
    pm <- c()
    if (L >= 1) pm <- c(pm, pout$p1)
    if (L >= 2) pm <- c(pm, pout$p2)
    if (L >= 3) pm <- c(pm, pout$p3)
    if (L >= 4) pm <- c(pm, pout$p4)
    if (L >= 5) pm <- c(pm, pout$p5)
    po <- .stk.organize(pm, stk$pos, t, derive.mu, fix.convo, fix.beta.a, fix.nu0.ratio) # po is the oraganized output
    names(po) <- c("mu", "nu0", "theta", "convo", "beta.a")
    do <- qsld(t = t, mu = po[1], nu0 = po[2], theta = po[3], convo = po[4], beta.a = po[5])
    if (derive.mu) {
        do@mu <- po[1] <- unname(mean_x - sld.mean(do))
    }
    if (verbose.interval > 0) {
        print("output from optimx:")
        print(list(pm=pm, po=po))
    }
    

    return(list(qsld=do, optimx.out=pout))
}
### <---------------------------------------------------------------------->
.stk.push <- function(stk, p, v) {
    stk$pm <- c(stk$pm, v)
    if (is.na(stk$pos[p])) stk$pos[p] <- length(stk$pm)
    else stop(paste("push failed, pos already occupied:", p))
    stk
}
.stk.from_qsld <- function(stk, d, derive.mu, fix.convo, fix.beta.a, fix.nu0.ratio) {
    if (! derive.mu) stk <- .stk.push(stk, 1, d@mu/d@t)
    if (is.na(fix.nu0.ratio)) stk <- .stk.push(stk, 2, d@nu0*100)
    stk <- .stk.push(stk, 3, d@theta*100)
    if (is.na(fix.convo)) stk <- .stk.push(stk, 4, d@convo)
    if (is.na(fix.beta.a)) stk <- .stk.push(stk, 5, d@beta.a)
    return(stk)
}
.stk.pop <- function(stk, p, factor=1) {
    if (is.na(stk$pos[p])) stop(paste("pop failed, pos without pointer:", p))
    L <- length(stk$pm)
    if (L==0) stop(paste("pop failed, pm has no more data"))
    stk$po[p] <- stk$pm[1]*factor
    stk$pos[p] <- NaN
    stk$pm <- if (L >= 2) stk$pm[2:L] else numeric(0)
    stk
}
.stk.organize <- function(pm, pos, t, derive.mu, fix.convo, fix.beta.a, fix.nu0.ratio) {
    # internal pop stack structure
    s <- list( 
        pm = pm # pm is raw input parameter array
        ,po = rep(NaN, 5) # po is the organized parameter array
        ,pos = pos # position of parameters, in the order of mu, nu0, theta, convo, beta.a
    )
    if (derive.mu) s$po[1] <- 0 # place holder for mu
    else s <- .stk.pop(s, 1, t) # mu
    
    if (is.na(fix.nu0.ratio)) {
        s <- .stk.pop(s, 2, 1/100) # nu0
        s <- .stk.pop(s, 3, 1/100) # theta
    }
    else {
        s <- .stk.pop(s, 3, 1/100) # theta
        s$po[2] <- s$po[3] * fix.nu0.ratio # nu0 = theta * fix.nu0.ratio
    }
    if (is.na(fix.convo)) s <- .stk.pop(s, 4) # convo
    else s$po[4] <- fix.convo 
    
    if (is.na(fix.beta.a)) s <- .stk.pop(s, 5) # beta.a
    else s$po[5] <- fix.beta.a
    
    if (length(s$pm) > 0) {
        msg <- "organize utilility failed, parameter array is not consumed entirely:"
        stop(sprintf("%s %s", msg, s))
    }
    return(s$po)
}
### <---------------------------------------------------------------------->

