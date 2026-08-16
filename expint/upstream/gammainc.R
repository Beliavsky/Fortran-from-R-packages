### == expint: Exponential Integral and Incomplete Gamma Function ==
###
### The incomplete gamma function
###
###    G(a,x) = int_x^infty t^{a-1} exp(-t) dt
###
### for a *real* and x >= 0. Note the order of the arguments.
###
### AUTHOR: Vincent Goulet <vincent.goulet@act.ulaval.ca>

gammainc <- function(a, x)
    .External(C_expint_do_gammainc, a, x)
