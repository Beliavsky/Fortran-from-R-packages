## Not exported, used to make  'R CMD check <pkg>'  be faster *or* more extensive:
doExtras <- function(int = interactive()) {
    int || nzchar(Sys.getenv("R_DPQ_check_extra")) ||
        identical("true", unname(Sys.getenv("R_PKG_CHECKING_doExtras")))
}

## Convenient argument checking

all_mpfr <- function(...) {
    ## Remain lazy ==> do *NOT* use list(...)  which evaluates all
    for(i in seq_len(...length())) if(!inherits(...elt(i), "mpfr")) return(FALSE)
    ## else return
    TRUE
}

any_mpfr <- function(...) {
    ## Remain lazy ==> do *NOT* use list(...)  which evaluates all
    for(i in seq_len(...length())) if(inherits(...elt(i), "mpfr")) return(TRUE)
    ## else return
    FALSE
}

## MM: From ~/R/D/r-devel/R/src/nmath/pgamma.c : --- Help page now  >>> ../man/logspace.add.Rd <<<<
logspace.add <- function(lx, ly) pmax(lx, ly) + log1p(exp(-abs(lx - ly)))

logspace.sub <- function(lx, ly) lx + log1mexp(lx - ly)

## Calling C in ../src/DPQ-misc.c :
log1mexpC <- function(x) .Call(C_R_log1mexp, x)
log1pexpC <- function(x) .Call(C_R_log1pexp, x)
log1pmxC  <- function(x) .Call(C_R_log1pmx,  x) # < TODO: add tuning args, as in log1pmx()
lgamma1pC <- function(x) .Call(C_R_lgamma1p, x)

pow <- function(x, y, try.int.y = TRUE) {
    if(try.int.y && ((Iy <- is.integer(y)) ||
                     all(abs(y) <= .Machine$integer.max)))
        pow_di(x, if(Iy) y else as.integer(y))
    else
        .Call(C_dpq_pow, x, y)
}
## for low-level comparisons (e.g. speed):
.pow   <- function(x, y) .Call(C_dpq_pow, x, y)
pow_di <- function(x, y) .Call(C_dpq_pow_di, x, y)

lgammacor <- function (x, nalgm = 5, xbig = 2^26.5)
    ## Hardwired in R's C code:  nalgm = 5, xbig = 2^26.5 = 94906265.62425156
    .Call(C_R_lgammacor, x, nalgm, xbig)

### was as 'form01.prec' in  source("~/R/MM/MISC/Util.R")
format01prec <- function(x, digits = getOption("digits"), width = digits + 2,
                        eps = 1e-6, ...,
                        FUN = function(x,...) formatC(x, flag='-',...))
{
  ## Purpose: format numbers in [0,1] with "precise" result,
  ##          using "1-.." if necessary.
  ## -------------------------------------------------------------------------
  ## Arguments: x:      numbers in [0,1]; (still works if not)
  ##            digits, width: number of digits, width to use with 'FUN'
  ##            eps:    Use '1-' iff  x in  (1-eps, 1] -- 1e-6 is OPTIMAL
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 14 May 1997
  if(as.integer(digits) < 4) stop('digits must be >= 4')
  if(eps < 0 || eps > .1) stop('eps must be in [0, .1]')
  i.swap <- 1-eps < x  &  x <= 1 #-- Use "1- ." if <~ 1,  normal 'FUN' otherwise
  r <- character(length(x))
  if(hasNA <- anyNA(i.swap)) {
    iNA <- is.na(i.swap)
    i.swap  <-  i.swap & !iNA
    ni.swap <- !i.swap & !iNA
  } else
    ni.swap <- !i.swap
  if(any(i.swap))
      r[i.swap] <- paste0("1-",
                          FUN(1-x[i.swap],
                              ## -5: '1-' + 4 for exponent -1 for '0' (in other case)
                              digits=digits - 5, width=width-2, ...))
  if(any(ni.swap))
    r[ni.swap] <- FUN(x[ni.swap], digits=digits, width=width,...)
  if(hasNA)
    r[iNA] <- "NA"
  attributes(r) <- attributes(x)
  r
}

## More calling C in ../src/DPQ-misc.c :

## C:  r = frexp (x, &e); // => r in  [0.5, 1) and 'e' (int) such that  x = r * 2^e
frexp <- function(x) .Call(C_R_frexp, x)

ldexp <- function(f, E) .Call(C_R_ldexp, f, E) # // ldexp(f, E) := f * 2^E

modf <- function(x) .Call(C_R_modf, x)

### Chebyshev Polynomial Evaluation ------- ../src/chebyshev.c ---------
### interface to R's Mathlib chebychev_init() and chebychev_eval()

##' Returns n_a such that coef[1+(0:n_a)] should be used[] (for max.error eta)
chebyshev_nc <- function(coef, eta = .Machine$double.eps/20) {
    stopifnot(length(eta) == 1, eta >= 0, length(coef) >= 1)
    .Call(C_R_chebyshev_nt, coef, eta)
}

##' Returns the \sum_{j=0}^n c_j T_j(x)  where c_0 = 2*coef[1], c_j := coef[j+1], j >= 1,
##'    n := nterms  and  T_j() is the Chebyshev polynomial of degree j
chebyshevEval <- function(x, coef,
                          nc = chebyshev_nc(coef, eta),
                          eta = .Machine$double.eps/20) {
    stopifnot(nc == (n <- as.integer(nc)), length(n) == 1L, n >= -1L)
    .Call(C_R_chebyshev_eval, x, coef, n)
}

.chebyshevEval <- function(x, coef, nc) .Call(C_R_chebyshev_eval, x, coef, nc)


##' Returns an \R  function(x)  which
chebyshevPoly <- function(coef,
                          nc = chebyshev_nc(coef, eta),
                          eta = .Machine$double.eps/20) {
    stopifnot(nc == (n <- as.integer(nc)), length(n) == 1L, n >= -1L)
    rm(nc, eta)
    function(x) .chebyshevEval(x, coef, n)
}
