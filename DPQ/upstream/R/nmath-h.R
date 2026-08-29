ISNAN <- is.na
ISNA <- function(x) is.na(x) && !is.nan(x)

ML_WARN_return_NAN <- function() { ML_WARNING(ME_DOMAIN, ""); return(ML_NAN); }

ML_VALID <- function(x)	!is.na(x) # C: !ISNAN(x)

ME_NONE      <- 0 # /*	no error */
ME_DOMAIN    <- 1 # /*	argument out of domain */
ME_RANGE     <- 2 # /*	value out of range */
ME_NOCONV    <- 4 # /*	process did not converge */
ME_PRECISION <- 8 # /*	does not have "full" precision */
ME_UNDERFLOW <- 16 #/*	and underflow occured (important for IEEE)*/

MATHLIB_WARNING <- function(msg, s) warning(warningCondition(msg, s,
							     class = "mathlib_warning",
							     call = sys.call(-2L)))
## /* For a long time prior to R 2.3.0 ML_WARNING did nothing.
##    We don't report ME_DOMAIN errors as the callers collect ML_NANs into
##    a single warning.
##  */
ML_WARNING <- function(x, s) {
    if(x > ME_DOMAIN) {
        msg = "";
        switch(x,
               ME_DOMAIN =
                   msg <- gettext("argument out of domain in '%s'n"),
               ME_RANGE =
                   msg <- gettext("value out of range in '%s'n"),
               ME_NOCONV =
                   msg <- gettext("convergence failed in '%s'n"),
               ME_PRECISION =
                   msg <- gettext("full precision may not have been achieved in '%s'n"),
               ME_UNDERFLOW =
                   msg <- gettext("underflow occurred in '%s'n"))
    }
    ## MATHLIB_WARNING(msg, s); :
    warning(warningCondition(sprintf(msg, s),
			     call = sys.call(-1L), ME_kind = x, class = "ML_WARNING"))
}

## C: # define MATHLIB_WARNING(fmt,x)		warning(fmt,x)
