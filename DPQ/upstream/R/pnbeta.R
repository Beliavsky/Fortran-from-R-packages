pnbetaAppr2v1 <- ## pnbeta.appr2 <-
    function(x, a, b, ncp = 0, lower.tail=TRUE, log.p=FALSE)
{
  ## Purpose: "Approximation 2" of
  ## Chattamvelli, R. \& Shanmugam, F.  (1997), `Algorithm {AS 310}: ...
  ## JRSS C (Applied Statistics)}  46(1), 146--156; notably p.155
  ## ----------------------------------------------------------------------
  ## Arguments: as for pbeta(), but named more traditionally
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 16 Oct 2007, 11:36

                                        # λ <- ncp -- not liked by CRAN checks
    a2l <- 2*a + ncp
    c. <- (a2l/(2*a)) ^ (1/3)
    d <- ((b * x)/(a*(1-x))) ^ (1/3)
    ##        --      --- <==> allow to input "1-x" for x ~= 1
    ## sigma_1^2  and sigma_2 ^2 --> mu_j = 1 - sigma_j^2  (j = 1,2):
    s1.2 <- 4*(a + ncp) / (9 * a2l^2)
    s2.2 <- 1/(9 * b)
    mu.L <- c. * (1 - s1.2) - d   * (1 - s2.2)
    si.L <- sqrt(c.^2 * s1.2     + d^2 * s2.2)
    pnorm( - mu.L / si.L, lower.tail=lower.tail, log.p = log.p)
}

pnbetaAppr2 <- function(x, a, b, ncp = 0, lower.tail=TRUE, log.p=FALSE)
{
  ## Purpose: "Approximation 2" of
  ## Chattamvelli, R. \& Shanmugam, F.  (1997), `Algorithm {AS 310}: ...
  ## JRSS C (Applied Statistics)}  46(1), 146--156; notably p.155
  ## ----------------------------------------------------------------------
  ## Arguments: as for pbeta(), but named more traditionally
  ## ----------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 16 Aug 2018;
  ##         *Math* simplified: c == 1 (+ cancellation of "/a")
                                        # λ <- ncp
    a2l <- 2*a + ncp
    d <- ((2*b * x)/(a2l*(1-x))) ^ (1/3)
    ##        --      --- <==> allow to input "1-x" for x ~= 1
    ## sigma_1^2  and sigma_2 ^2 --> mu_j = 1 - sigma_j^2  (j = 1,2):
    s1.2 <- 4*(a + ncp) / (9 * a2l^2)
    s2.2 <- 1/(9 * b)
    mu.L <- (1 - s1.2) - d * (1 - s2.2)
    si.L <- sqrt(s1.2  + d^2 * s2.2)
    pnorm( - mu.L / si.L, lower.tail=lower.tail, log.p = log.p)
}

pnbetaAS310 <- function(x, a, b, ncp = 0, lower.tail=TRUE, log.p=FALSE,
                        useAS226 = (ncp < 54.), errmax = 1e-6, itrmax = 100)
{
    stopifnot(length(lower.tail) == 1, length(log.p) == 1,
              length(errmax) == 1, length(itrmax) == 1,
              length(useAS226) %in% c(1L, length(ncp)))
    ## Two cases:
    ##  1)  length(a) == length(b) == length(ncp) == n := length(x)   <<  all_n <- TRUE
    ##  2)  length(a) == length(b) == length(ncp) == 1                <<  all_n <- FALSE
    n <- length(x <- as.numeric(x))
    all_n  <- length(a) > 1
    if(all_n) {
        ## recycle the first 4 arguments to common length
        ## [... yes, would be nicer with .Call() where the C code could nicely "wrap around" via (i % n) -- FIXME?
        n <- max(n,
                 length(a <- as.numeric(a)),
                 length(b <- as.numeric(b)),
                 length(ncp <- as.numeric(ncp)))
        if(n != length(x)) x <- rep_len(x, n)
        if(n != length(a)) a <- rep_len(a, n)
        if(n != length(b)) b <- rep_len(b, n)
        if(n != length(ncp)) { ncp <- rep_len(ncp, n); useAS226 <- rep_len(useAS226, n) }
    } else { ## only x[] is of length n .. the other three of length 1
        stopifnot(length(a <- as.numeric(a)) == 1,
                  length(b <- as.numeric(b)) == 1,
                  length(ncp <- as.numeric(ncp)) == 1)
    }
    r <- .C(C_ncbeta, # ../src/310-pnbeta.c
            a, b, ncp, x,
            as.integer(n),
            as.logical(useAS226),
            as.double(errmax),
            as.integer(itrmax),
            ifault = as.integer(all_n),# input/output
            res = double(n))[c("ifault","res")]
    if(r$ifault) ## TODO: switch(r$ifault, ....) for different error messages
        stop(sprintf("ifault=%d from error in C code ncbeta()", r$ifault))
    else
        r$res
}
