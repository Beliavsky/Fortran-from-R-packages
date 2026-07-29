## By Marius Hofert

### Composite distributions -- examples ########################################

library(qrmtools)

## For now:
dcomposite <- qrmtools:::dcomposite
pcomposite <- qrmtools:::pcomposite
qcomposite <- qrmtools:::qcomposite
rcomposite <- qrmtools:::rcomposite


### Example 1: (N(2,2), Exp(3))

## Setup
x <- seq(-5, 5, length.out = 256)
cuts <- 3
weights <- rep(1/2, 2)
distr <- list(list("norm", c(2,2)), list("exp", 3))
## d*()
dc <- dcomposite(x, cuts = cuts, distr = distr, weights = weights)
plot(x, dc, type = "l", ylab = "dcomposite(x)")
## p*()
pc <- pcomposite(x, cuts = cuts, distr = distr, weights = weights)
plot(x, pc, type = "l", ylab = "pcomposite(x)")
## q*()
p <- seq(0, 1, length.out = 256)
qc <- qcomposite(p, cuts = cuts, distr = distr, weights = weights)
plot(p, qc, type = "l", ylab = "qcomposite(p)")
## r*()
set.seed(271)
rc <- rcomposite(1000, cuts = cuts, distr = distr, weights = weights)
plot(rc)
plot(pcomposite(rc, cuts = cuts, distr = distr, weights = weights))
## should be U[0,1] => looks good


### Example 2: (N(2,2), Par(3))

## Setup
x <- seq(-5, 5, length.out = 256)
cuts <- 3
weights <- rep(1/2, 2)
dPar <- function(x, theta = 3) theta*(1+x)^(-theta-1)
pPar <- function(x, theta = 3) 1-1/(1+x)^theta
qPar <- function(x, theta = 3) (1-x)^(-1/theta)-1
## d*()
distr <- list(list("norm", c(2,2)), list(pPar, dPar))
dc <- dcomposite(x, cuts = cuts, distr = distr, weights = weights)
plot(x, dc, type = "l", ylab = "dcomposite(x)")
## p*()
## Note: it's fine here to use the same 'distr' (distr[[2]][[2]] is not used here)
pc <- pcomposite(x, cuts = cuts, distr = distr, weights = weights)
plot(x, pc, type = "l", ylab = "pcomposite(x)")
## q*()
## Note: now we have to provide the quantile function
distr <- list(list("norm", c(2,2)), list(pPar, qPar))
p <- seq(0, 1, length.out = 256)
qc <- qcomposite(p, cuts = cuts, distr = distr, weights = weights)
plot(p, qc, type = "l", ylab = "qcomposite(p)")
## r*()
set.seed(271)
rc <- rcomposite(1000, cuts = cuts, distr = distr, weights = weights)
plot(rc)
plot(pcomposite(rc, cuts = cuts, distr = distr, weights = weights))
## Should be U[0,1]  => looks good


### Example 3: (ECDF, Par(3))

## Setup
set.seed(271)
Y <- rgamma(100, shape = 1/2)
## d*()
distr <- list(ecdf(Y), list(pPar, dPar))
dc <- dcomposite(x, cuts = cuts, distr = distr, weights = weights)
plot(x, dc, type = "l", ylab = "dcomposite(x)") # note: we don't have a density on [0,3] here
## p*()
## Note: it's fine here to use the same 'distr' (distr[[2]][[2]] is not used here)
pc <- pcomposite(x, cuts = cuts, distr = distr, weights = weights)
plot(x, pc, type = "l", ylab = "pcomposite(x)")
## q*()
## Note: now we have to provide the (empirical) quantile function
distr <- list(list(ecdf(Y), function(p) quantile(x = Y, probs = p)), list(pPar, qPar))
p <- seq(0, 1, length.out = 256)
qc <- qcomposite(p, cuts = cuts, distr = distr, weights = weights)
plot(p, qc, type = "l", ylab = "qcomposite(p)")
## r*()
set.seed(271)
rc <- rcomposite(1000, cuts = cuts, distr = distr, weights = weights)
plot(rc)
plot(pcomposite(rc, cuts = cuts, distr = distr, weights = weights))
## Should be U[0,1] => looks good
