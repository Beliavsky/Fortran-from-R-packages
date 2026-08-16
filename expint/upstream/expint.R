### == expint: Exponential Integral and Incomplete Gamma Function ==
###
### Functions to compute the exponential integral Ei and the
### exponential integral of order 'n = 1, 2, ...'
###
###     E_n = int_x^infty exp(-t)/t^n dt.
###
### Function 'expint' is a unified and fully vectorized interface with
### default the most common case E_1. The other functions are simpler,
### slightly faster interfaces to E_1, E_2, E_n and Ei.
###
### When 'scale' is TRUE, the value returned is scaled by exp(x).
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

expint <- function(x, order = 1L, scale = FALSE)
    .External(C_expint_do_expint, "En", x, order, scale)

expint_E1 <- function(x, scale = FALSE)
    .External(C_expint_do_expint, "E1", x, scale)

expint_E2 <- function(x, scale = FALSE)
    .External(C_expint_do_expint, "E2", x, scale)

expint_En <- function(x, order, scale = FALSE)
    .External(C_expint_do_expint, "En", x, order[1L], scale)

expint_Ei <- function(x, scale = FALSE)
    -.External(C_expint_do_expint, "E1", -x, scale)
