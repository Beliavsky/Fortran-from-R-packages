
library("mvtnorm")

options(digits = 5)
tol <- sqrt(sqrt(.Machine$double.eps))

set.seed(29078)

EVAL <- function(...) {}

if (require("numDeriv", quietly = TRUE))
    EVAL <- eval

chk <- function(...) stopifnot(isTRUE(all.equal(..., check.attributes = FALSE, 
                                                tol = tol)))

x <- mvnorm()
J <- 3
M <- diag(1:J)
rownames(M) <- colnames(M) <- LETTERS[1:J]
(x <- mvnorm(mean = runif(J), chol = M))
margDist(x, which = 2:J)
margDist(x, which = 2:J)
condDist(x, which_given = 1, given = matrix(1))

logLik(x, obs = M[-1,-1])
logLik(margDist(x, which = 2:J), obs = M[-1,-1])


thischeck <- expression({

#set.seed(29)

l <- matrix(pl <- runif(J * (J - 1) / 2 * Ns),
            ncol = Ns)
colnames(l) <- paste0("i", 1:Ns)
L <- ltMatrices(l, byrow = BYROW, names = LETTERS[1:J])
m <- matrix(pm <- rnorm(J * Nm), ncol = Nm)
rownames(m) <- LETTERS[1:J]
colnames(m) <- paste0("i", 1:Nm)
if (CHOL) {
    x <- mvnorm(mean = m, chol = L)
} else {
    x <- mvnorm(mean = m, invchol = L)
}
obs <- simulate(x, nsim = N, standardize = TRUE)

tfun <- function(parm, perm = LETTERS[1:J], FUN = "logLik", chol = TRUE, ...) {
    args <- list(...)
    p1 <- parm[1:length(pm)]
    p2 <- parm[length(pm) + 1:length(pl)]
    p3 <- parm[-(1:(length(pm) + length(pl)))] 
    p2 <- matrix(p2, ncol = Ns)
    L <- ltMatrices(p2, names = LETTERS[1:J])
    L <- ltMatrices(p2, byrow = BYROW, names = LETTERS[1:J])
    m <- matrix(p1, ncol = Nm)
    rownames(m) <- LETTERS[1:J]
    obs <- matrix(p3, ncol = N)
    rownames(obs) <- LETTERS[1:J]
    if (chol) {
        x <- mvnorm(mean = m, invchol = L)
    } else {
        x <- mvnorm(mean = m, chol = L)
    }
    args$object <- x
    args$obs <- obs[perm,,drop = FALSE]
    do.call(FUN, args)
}

ll <- tfun
sc <- function(...) tfun(..., FUN = "lLgrad")

l1 <- logLik(x, obs)
#lLgrad(x, obs)

prm <- c(pm, pl, obs)
ll(prm)
s <- sc(prm)
sa <- c(if (Nm > 1) s$mean else rowSums(s$mean),
        if (Ns > 1) Lower_tri(s$scale) else rowSums(Lower_tri(s$scale)), 
        s$obs)
sn <- grad(ll, prm)
chk(sa, sn)

l2 <- logLik(x, obs, standardize = TRUE)
# lLgrad(x, obs, standardize = TRUE)

ll(prm, standardize = TRUE)
s <- sc(prm, standardize = TRUE)
sa <- c(if (Nm > 1) s$mean else rowSums(s$mean),
        if (Ns > 1) Lower_tri(s$scale) else rowSums(Lower_tri(s$scale)), 
        s$obs)
sn <- grad(ll, prm, standardize = TRUE)
chk(sa, sn)

l1p <- logLik(x, obs = obs[perm <- sample(rownames(obs)),,drop = FALSE])
#lLgrad(x, obs = obs[perm,,drop = FALSE])

chk(l1, l1p)

ll(prm, perm = perm)
s <- sc(prm, perm = perm)
sa <- c(if (Nm > 1) s$mean else rowSums(s$mean),
        if (Ns > 1) Lower_tri(s$scale) else rowSums(Lower_tri(s$scale)), 
        s$obs[LETTERS[1:J],])
sn <- grad(ll, prm, perm = perm)
chk(sa, sn)

l2p <- logLik(x, obs = obs[perm,,drop = FALSE], standardize = TRUE)
# lLgrad(x, obs = obs[perm,,drop = FALSE], standardize = TRUE)

chk(l2, l2p)

ll(prm, perm = perm, standardize = TRUE)
s <- sc(prm, perm = perm, standardize = TRUE)
sa <- c(if (Nm > 1) s$mean else rowSums(s$mean),
        if (Ns > 1) Lower_tri(s$scale) else rowSums(Lower_tri(s$scale)), 
        s$obs[LETTERS[1:J],])
sn <- grad(ll, prm, perm = perm, standardize = TRUE)
chk(sa, sn)

logLik(x, obs = obs[perm[-1],,drop = FALSE])
# lLgrad(x, obs = obs[perm,,drop = FALSE])

ll(prm, perm = perm[-1])
s <- sc(prm, perm = perm[-1])
s$obs <- rbind(s$obs, 0)
rownames(s$obs)[nrow(s$obs)] <- perm[1]
sa <- c(if (Nm > 1) s$mean else rowSums(s$mean),
        if (Ns > 1) Lower_tri(s$scale) else rowSums(Lower_tri(s$scale)), 
        s$obs[LETTERS[1:J],])
sn <- grad(ll, prm, perm = perm[-1])

logLik(x, obs = obs[perm[-1],,drop = FALSE], standardize = TRUE)
# lLgrad(x, obs = obs[perm[-1],,drop = FALSE], standardize = TRUE)

ll(prm, perm = perm[-1], standardize = TRUE)
s <- sc(prm, perm = perm[-1], standardize = TRUE)
s$obs <- rbind(s$obs, 0)
rownames(s$obs)[nrow(s$obs)] <- perm[1]
sa <- c(if (Nm > 1) s$mean else rowSums(s$mean),
        if (Ns > 1) Lower_tri(s$scale) else rowSums(Lower_tri(s$scale)), 
        s$obs[LETTERS[1:J],])
sn <- grad(ll, prm, perm = perm[-1], standardize = TRUE)
chk(sa, sn)
})

J <- 4
Ns <- 1
Nm <- 1
N <- 3 #max(Ns, Nm)
BYROW <- FALSE
CHOL <- FALSE
EVAL(thischeck)

J <- 4
Ns <- 3
Nm <- 3
N <- max(Ns, Nm)
BYROW <- TRUE
CHOL <- TRUE
EVAL(thischeck)

J <- 4
Ns <- 1
Nm <- 3
N <- max(Ns, Nm)
BYROW <- TRUE
CHOL <- TRUE
EVAL(thischeck)

J <- 4
Ns <- 3
Nm <- 1
N <- max(Ns, Nm)
BYROW <- TRUE
CHOL <- TRUE
EVAL(thischeck)

J <- 4
Ns <- 3
Nm <- 1
N <- max(Ns, Nm)
BYROW <- TRUE
CHOL <- FALSE
EVAL(thischeck)

J <- 4
Ns <- 3
Nm <- 3
N <- max(Ns, Nm)
BYROW <- FALSE
CHOL <- TRUE
EVAL(thischeck)

J <- 4
Ns <- 3
Nm <- 3
N <- max(Ns, Nm)
BYROW <- FALSE
CHOL <- FALSE
EVAL(thischeck)

J <- 4
Ns <- 3
Nm <- 3
N <- max(Ns, Nm)
BYROW <- FALSE
CHOL <- FALSE
EVAL(thischeck)

