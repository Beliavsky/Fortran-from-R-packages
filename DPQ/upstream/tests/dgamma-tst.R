library(DPQ)

set.seed(17)
x <- rlnorm(500)
dg.R.1  <- sapply(x, dgamma.R, shape = 11.111)
dg.R.1L <- sapply(x, dgamma.R, shape = 11.111, log = TRUE)
## print observed differences
all.equal(dgamma  (x, 11.111          ), dg.R.1,  tolerance = 0) # 7.844e-16
all.equal(dgamma  (x, 11.111, log=TRUE), dg.R.1L, tolerance = 0) # 2.148e-16
stopifnot(exprs = {
    all.equal(dg.R.1,  dgamma.R(x, 11.111          ), tolerance = 0) # just vectorization
    all.equal(dg.R.1L, dgamma.R(x, 11.111, log=TRUE), tolerance = 0) #  "    "
    all.equal(dg.R.1,  dgamma(x, 11.111          ), tolerance = 1e-14)
    all.equal(dg.R.1L, dgamma(x, 11.111, log=TRUE), tolerance = 1e-14)
    all.equal(dgamma  (x, 14.99), tolerance = 2e-15,
              dgamma.R(x, 14.99, dpois_r_args = list(verbose=TRUE) -> dgR149))
})



