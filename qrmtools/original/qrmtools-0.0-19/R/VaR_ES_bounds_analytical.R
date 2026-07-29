### Tools for computing the VaR_alpha and ES_alpha bounds ######################

### 1 Crude VaR bounds (for both best and worst VaR) ###########################

##' @title Crude bounds for any VaR_alpha
##' @param level confidence level
##' @param qF (list of) marginal quantile functions
##' @param d if qF is a function, dimension for the homogeneous case; ignored
##'        otherwise
##' @param ... ellipsis argument passed to qF()
##' @return 2-vector containing crude VaR_alpha bounds
##' @author Marius Hofert
##' @note See Lemma 2.1 in "Improved Algorithms for Computing Worst VaR"
crude_VaR_bounds <- function(level, qF, d = NULL, ...)
{
    ## ... are passed to *all* qF()
    if(is.list(qF)) { # inhomogeneous case
        if(!all(unlist(lapply(qF, is.function))))
            stop("qF has to be a (quantile) function or list of such")
        d <- length(qF)
        qF.low <- sapply(qF, function(qF.) qF.(level/d, ...))
        qF.up  <- sapply(qF, function(qF.) qF.((d-1+level)/d, ...))
        d * c(min(qF.low), max(qF.up))
    } else { # homogeneous case
        if(!is.function(qF))
            stop("qF has to be a (quantile) function or list of such")
        if(is.null(d))
            stop("if qF is a single function, d has to be a positive integer")
        qF.low <- qF(level/d, ...)
        qF.up  <- qF((d-1+level)/d, ...)
        d * c(qF.low, qF.up)
    }
}


### 2 Explicit worst/best VaR in the homogeneous case ##########################

### Dual bound #################################################################

##' @title D(s,t) = d \int_{t}^{s-(d-1)t} \bar{F}(x) dx / (s-dt)
##' @param s real number
##' @param t real number < s/d
##' @param d dimension (integer > 2)
##' @param pF marginal distribution function (same for all d)
##' @param ... ellipsis argument passed to integrate()
##' @return D(s,t)
##' @author Marius Hofert
##' @note If t -> s/d-, l'Hospital's Rule shows that D(s, s/d) = d\bar{F}(s/d)
dual_bound_2 <- function(s, t, d, pF, ...)
{
    stopifnot(length(t) == 1)
    if(t > s/d) stop("t must be <= s/d")
    ## use D(s,t) = d( 1-\int_{t}^{s-(d-1)t} F(x) dx/(s-d*t) ) in this case
    if(t == s/d) d*(1-pF(s/d)) else
    d * (1 - (1/(s-d*t)) * integrate(pF, lower = t, upper = s-(d-1)*t, ...)$value)
}

##' @title Auxiliary function \bar{F}(t) + (d-1) * \bar{F}(s-(d-1)*t)
##' @param s real number
##' @param t real number < s/d
##' @param d dimension (integer > 2)
##' @param pF marginal distribution function (same for all d)
##' @return \bar{F}(t) + (d-1) * \bar{F}(s-(d-1)*t)
##' @author Marius Hofert
dual_bound_2_deriv_term <- function(s, t, d, pF)
    1-pF(t) + (d-1)*(1-pF(s-(d-1)*t))

##' @title Dual bound D(s)
##' @param s real number
##' @param d dimension (integer > 2)
##' @param pF marginal distribution function (same for all d)
##' @param ... ellipsis argument passed to dual_bound_2()'s integrate()
##' @return D(s)
##' @author Marius Hofert
##' @note The "first-order condition" (second equality in (14) in 2)) comes from the
##'       fact that
##'       (d/dt) D(s,t) = [ (-d)[\bar{F}(s-(d-1)t)(d-1)+\bar{F}(t)](s-dt) +
##'                         d^2 \int_{t}^{s-(d-1)t} \bar{F}(x) dx ] / (s-dt)^2 = 0
##'       if and only if
##'       d (\int_{t}^{s-(d-1)t} \bar{F}(x) dx) / (s-dt) = \bar{F}(s-(d-1)t)(d-1)-\bar{F}(t)
##'       => solving d (\int_{t}^{s-(d-1)t} \bar{F}(x) dx) / (s-dt) -
##'                  (\bar{F}(s-(d-1)t)(d-1)-\bar{F}(t)) = 0
##'          as a function in t for sufficiently large s leads to D(s)
dual_bound <- function(s, d, pF, tol = .Machine$double.eps^0.25, ...)
{
    stopifnot(length(s) == 1, s >= 0)
    if(s > 0) {
        ## h(s, t)
        h <- function(t) dual_bound_2(s, t = t, d = d, pF = pF, ...) -
            dual_bound_2_deriv_term(s, t = t, d = d, pF = pF)
        ## Note: h(t) -> 0 for t -> s/d- which is bad for uniroot() as the
        ##       latter will simply stop with the root t = s/d => we thus set f.upper > 0
        h.up <- -h(0) # guarantee that uniroot() doesn't fail due to root s/d
        t. <- uniroot(h, interval = c(0, s/d), f.upper = h.up, tol = tol)$root # optimal t in Equation (12) [= arginf]
        dual_bound_2_deriv_term(s, t = t., d = d, pF = pF) # dual bound D(s) in Equation (12) [= inf]
    } else {
        ## If s = 0, then t in [0, s/d] requires t to be 0 *and* f(0) = 0, so
        ## 0 is a root (as s/d). Furthermore, at t = 0 (and with s = 0),
        ## dual_bound_2_deriv_term(...) = d
        d
    }
}


### Wang's methods #############################################################

##' @title Scaled right-hand side term in the objective function for computing
##'        worst VaR as in McNeil, Frey, Embrechts (2015, Prop. 8.32)
##' @param c evaluation point
##' @param level confidence level alpha
##' @param d dimension d
##' @param method character string giving the method
##'        generic = numerical integration; Wang.Par = Pareto distibution
##' @param ... ellipsis argument containing shape (for method = "Wang.Par")
##'        or qF (for method = "generic")
##' @return Right-hand side term in Prop. 3.1
##' @author Marius Hofert
##' @note for the correct 'c', this is the conditional expectation
Wang_h_aux <- function(c, level, d, method = c("generic", "Wang.Par"), ...)
{
    ddd <- list(...)
    method <- match.arg(method)
    switch(method,
    "generic" = {
        qF <- ddd$qF # needs 'qF()'
        a <- level + (d-1)*c
        b <- 1-c
        qF(a)*(d-1)/d + qF(b)/d
    },
    "Wang.Par" = {
        ## We don't use qF(a)*(d-1)/d + qF(b)/d for qF(x) = qPar(x, shape = shape)
        ## here as qF(b) = qF(1-c) and 1-c == 1 for small c > 0 => numerically,
        ## qF(b) = Inf then.
        th <- ddd$shape # needs 'shape'
        t1 <- (1-level)/c-(d-1)
        (c^(-1/th)/d) * ((d-1)*t1^(-1/th) + 1) - 1 # checked (= qF(a)*(d-1)/d + qF(b)/d)
    },
    stop("Wrong method"))
}

##' @title Objective function for computing the worst VaR as in
##'        McNeil, Frey, Embrechts (2015, Prop. 8.32)
##' @param c evaluation point
##' @param level confidence level alpha
##' @param d dimension d
##' @param method character string giving the method
##' @param ... ellipsis argument passed to Wang_h_aux() and integrate()
##' @return objective function for computing the worst VaR
##' @author Marius Hofert
Wang_h <- function(c, level, d, method = c("generic", "Wang.Par"), ...)
{
    stopifnot(0 <= c, c <= (1-level)/d) # sanity check (otherwise b > a)
    method <- match.arg(method)
    ddd <- list(...)

    ## Compute \bar{I}(a, b) = 1/(b-a)\int_a^b qF(y) dy = (subs) IE[L|L\in [qF(a), aF(b)]]
    Ibar <- switch(method,
    "generic" = {
        qF <- ddd$qF # needs 'qF()'
        if(c == (1-level)/d) { # Properly deal with limit c = (1-alpha)/d
            qF(1-(1-level)/d)
        } else {
            a <- level + (d-1)*c
            b <- 1-c
            ddd$qF <- NULL # remove from '...'
            int <- function(...)
                integrate(qF, lower = a, upper = b, ...)$value / (b-a)
            do.call(int, ddd)
        }
    },
    "Wang.Par" = {
        th <- ddd$shape # needs 'shape'
        if(c == (1-level)/d) { # Properly deal with limit c = (1-alpha)/d
            ((1-level)/d)^(-1/th) - 1
        } else {
            t1 <- (1-level)/c-(d-1)
            t2 <- 1-level-d*c
            if(th == 1) log(t1)/t2 - 1
            else (th/(1-th))*c^(1-1/th)*(1-t1^(1-1/th))/t2 - 1
        }
    },
    stop("Wrong method"))

    ## Return
    Ibar - Wang_h_aux(c, level = level, d = d, method = method, ...)
}


### Main wrapper function for computing the best/worst VaR in the homogeneous case

## Assumptions:
## - d = 2: ultimately decreasing density (for x >= x0), alpha >= F(x0)
## - "Wang": F needs to live on [0, Inf), admitting a positive density which is
##           ultimately decreasing (for x >= x0), alpha >= F(x0)
## - "dual": F needs to be continuous with unbounded support and and ultimately
##           decreasing density, F(0) = 0 (otherwise, 0 as a lower bound for
##           uniroot() in dual_bound() is not valid)

##' @title Compute the best/worst VaR_\alpha in the homogeneous case with:
##'        1) d = 2: Embrechts, Puccetti, Rueschendorf (2013, Proposition 2)
##'        2) d >= 3:
##'           "Wang": McNeil, Frey, Embrechts (2015, Prop. 8.32)
##'                   Integral evaluated numerically; needs smaller default
##'                   tolerance for uniroot()!
##'           "Wang.Par": The same, just with explicit formula for the integral
##'                       in the Pareto case. Note that this requires to
##'                       extend the (theoretically correct) initial interval and
##'                       a smaller tolerance (see vignette).
##'           "dual": Embrechts, Puccetti, Rueschendorf (2013, Proposition 4)
##'                   Numerically less stable; no formula for best VaR known (=> NA)
##' @param level confidence level
##' @param d dimension
##' @param method character string giving the method
##' @param interval initial interval
##' @param tol uniroot() x-tolerance
##' @param ... ellipsis arguments passed to Wang_h()
##' @return (best VaR, worst VaR) in the homogeneous case
##' @author Marius Hofert
##' @note (*) Typos:
##'       - Wang, Peng, Yang (2013): best VaR wrong
##'       - Published version of Embrechts, Puccetti, Rueschendorf, Wang, Beleraj (2014, Prop. 3.1):
##'         Best VaR and worst VaR formulas wrong
##'       - Updated version of Embrechts, Puccetti, Rueschendorf, Wang, Beleraj (2014)
##'         on Ruodu's website: Correct worst VaR but still wrong best VaR (Eq. (3.4)).
##'       - Both best and worst VaR are correct in McNeil, Frey, Embrechts (2015, Prop. 8.32)
VaR_bounds_hom <- function(level, d, method = c("Wang", "Wang.Par", "dual"),
                           interval = NULL, tol = NULL, ...)
{
    stopifnot(0 < level, level < 1, d >= 2)
    method <- match.arg(method)

    ## Deal with d == 2 first ##################################################

    if(d == 2) { # See Embrechts, Puccetti, Rueschendorf (2013, Prop. 2)
        if(method == "Wang.Par") {
            shape <- NULL # make CRAN check happy
            if(!hasArg(shape))
                stop("The Pareto case requires the parameter 'shape'")
            th <- list(...)$shape
            stopifnot(length(th) == 1, th > 0) # check shape here
            qF <- function(p) qPar(p, shape = th)
            return( c((1-level)^(-1/th)-1, 2*(((1-level)/2)^(-1/th)-1)) )
        } else {
            qF <- NULL # make CRAN check happy
            if(!hasArg(qF))
                stop("The quantile function qF of F is required")
            qF <- list(...)$qF
            return(c(qF(level), 2*qF((1+level)/2)))
        }
    }

    ## Best VaR for d >= 3 #####################################################

    best <-
        switch(method,
               "Wang" = {
                   qF <- NULL # make CRAN check happy
                   if(!hasArg(qF))
                       stop("Method 'Wang' requires the quantile function qF of F")
                   ddd <- list(...)
                   qF <- ddd$qF # get qF()
                   ddd$qF <- NULL # remove from '...'
                   int <- function(...)
                       integrate(qF, lower = 0, upper = level, ...)$value / level
                   max((d-1)*qF(0)+qF(level), # See (*) above
                       d * do.call(int, ddd))
               },
               "Wang.Par" = {
                   shape <- NULL # make CRAN check happy
                   if(!hasArg(shape))
                       stop("Method 'Wang.Par' requires the parameter 'shape'")
                   th <- list(...)$shape
                   stopifnot(length(th) == 1, th > 0) # check shape here
                   Ibar <- if(th == 1) {
                       -log1p(-level) - level
                   } else {
                       ((1-level)^(1-1/th)-1)/(1-1/th) - level
                   }
                   max((d-1)*0 + (1-level)^(-1/th)-1, # See (*) above
                       d * Ibar)
               },
               "dual" = { # "dual" only provides worst VaR
                   NA
               },
               stop("Wrong method"))

    ## Worst VaR for d >= 3  ###################################################

    if(is.null(tol)) # use smaller tol
        tol <- if(method == "Wang" || method == "Wang.Par") 2.2204e-16 # MATLAB default
               else .Machine$double.eps^0.25 # uniroot() default
    worst <- switch(method,
           "Wang" = {

               ## Check qF()
               qF <- NULL # make CRAN check happy
               if(!hasArg(qF))
                   stop("Method 'Wang' requires the quantile function qF of F")
               ## Check 'interval'
               if(is.null(interval)) interval <- c(0, (1-level)/d)
               else {
                   if(interval[1] < 0) stop("interval[1] needs to be >= 0")
                   if(interval[1] > (1-level)/d) stop("interval[2] needs to be <= (1-level)/d")
                   if(interval[1] >= interval[2]) stop("interval[1] needs to be smaller than interval[2]")
               }

               ## Compute (and adjust) function values at endpoints
               h.low <- Wang_h(interval[1], level = level, d = d, ...)
               if(is.na(h.low))
                   stop("Objective function at interval[1] is NA or NaN. Provide a larger interval[1].")
               h.up <- -h.low # avoid that uniroot() fails due to 0 at upper interval endpoint

               ## Root-finding on 'interval'
               c. <- uniroot(function(c) Wang_h(c, level = level, d = d, ...),
                             interval = interval, f.lower = h.low, f.upper = h.up, tol = tol)$root
               d * Wang_h_aux(c., level = level, d = d, ...)

           },
           "Wang.Par" = {

               ## Critical here (see vignette): smaller tolerance and extending the initial interval

               ## Check 'shape'
               shape <- NULL # make CRAN check happy
               if(!hasArg(shape))
                   stop("Method 'Wang.Par' requires the parameter 'shape'")
               th <- list(...)$shape
               stopifnot(length(th) == 1, th > 0) # check shape here

               ## Compute uniroot() initial interval
               if(is.null(interval)) {
                   low <- if(th > 1) {
                       r <- (1-level)/((d/(th-1)+1)^th + d-1)
                       r/2 # adjustment to guarantee numerically that h is of opposite sign (required for very large shape)
                   } else if(th == 1) {
                       e <- exp(1)
                       (1-level)/((d+1)^(e/(e-1))+d-1)
                   } else {
                       r <- (1-th)*(1-level)/d
                       r/2 # adjustment to guarantee numerically that h is of opposite sign (required for very small shape)
                   }
                   up <- if(th == 1) {
                       (1-level)/(3*d/2-1)
                   } else {
                       (1-level)*(d-1+th)/((d-1)*(2*th+d))
                   }
                   interval <- c(low, up)
               } else {
                   if(interval[1] < 0) stop("interval[1] needs to be >= 0")
                   if(interval[1] > (1-level)/d) stop("interval[2] needs to be <= (1-level)/d")
                   if(interval[1] >= interval[2]) stop("interval[1] needs to be smaller than interval[2]")
               }

               ## Root-finding on 'interval'
               h <- function(c) Wang_h(c, level = level, d = d, method = "Wang.Par", ...)
               c <- uniroot(h, interval = interval, tol = tol)$root
               d * Wang_h_aux(c, level = level, d = d, method = "Wang.Par", shape = th)

           },
           "dual" = {

               pF <- NULL # make CRAN check happy
               if(!hasArg(pF))
                   stop("Method 'dual' requires the distribution function pF")
               if(!hasArg(interval))
                   stop("Method 'dual' requires an initial interval c(s_l, s_u) to be given")
               uniroot(function(s) dual_bound(s, d = d, tol = tol, ...) - (1-level),
                       interval = interval, tol = tol)$root # s interval
               ## Note: We can't pass arguments to the inner root-finding

           },
           stop("Wrong method"))

           ## Return
           c(best, worst)
}
