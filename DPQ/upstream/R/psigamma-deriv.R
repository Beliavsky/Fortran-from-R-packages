
## SEXP R_dpsifn(SEXP x_, SEXP m_, SEXP deriv_1, SEXP kode2_)
##  result: a list(dpsi = <m x n> matrix of dpsi() values,
##                 nz   =  n-vector,
##                 ierr =  n-vector)
dpsifn <- function(x, m, deriv1 = 0L, k2 = FALSE) {
    stopifnot(is.numeric(x),
              m == as.integer(m), length(m) == 1, m >= 1,
              deriv1 == (d1 <- as.integer(deriv1)), length(d1) == 1, 0 <= d1, d1 <= 100,
              is.logical(k2), length(k2) == 1, !is.na(k2))
    r <- .Call(C_R_dpsifn, x, as.integer(m), as.integer(deriv1), k2)
    structure(r$dpsi
            , underflow = if(any(r$nz  )) r$nz   else FALSE
            , errorCode = if(any(r$ierr)) r$ierr else FALSE
              )
}
