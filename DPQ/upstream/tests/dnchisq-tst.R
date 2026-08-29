library(DPQ)

set.seed(1)
for(x in rlnorm(500))
    stopifnot(all.equal(dnchisqR(x,10, 2),
                        dchisq  (x,10, 2), tol=1e-7))

