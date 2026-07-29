
library("mvtnorm")

options(digits = 3)
tol <- 1e-1

set.seed(29)

EVAL <- function(...) {}

if (require("numDeriv", quietly = TRUE))
    EVAL <- eval

chk <- function(...) stopifnot(isTRUE(all.equal(..., check.attributes = FALSE, tol = sqrt(sqrt(.Machine$double.eps)))))

thischeck <- expression({
J <- 5
p <- sample(1:J)
if (isTRUE(all.equal(p, 1:J)))
    warning("Checks for id permutation meaningless")
P <- matrix(0, nrow = J, ncol = J)
P[cbind(1:J, p)] <- 1

L <- as.invchol(ltMatrices(1 + runif(J * (J + 1) / 2), diag = TRUE, byrow = BYROW))
mL <- as.array(L)[,,1]
S <- invchol2cov(L)
mS <- as.array(S)[,,1]
mSp <- mS[p,p]

chk(P %*% mS %*% t(P), mSp)

O <- invchol2pre(L)
mO <- as.array(O)[,,1]
chk(solve(P %*% mO %*% t(P)), mSp)
chk(solve(P %*% t(mL) %*% mL %*% t(P)), mSp)

C <- invchol2chol(L)
mC <- as.array(C)[,,1]
chk(P %*% mC %*% t(mC) %*% t(P), mSp)

Ct <- t(chol(mS[p,p]))
chk(Ct %*% t(Ct), mSp)

chk(as.array(invchol2cov(aperm(L, perm = p)))[,,1], mSp)
chk(as.array(chol2cov(aperm(C, perm = p)))[,,1], mSp)

N <- 10000
obs <- matrix(rnorm(J * N), ncol = N)
obs <- Mult(C, obs)

ll1 <- ldmvnorm(obs = obs, chol = C)
ll2 <- ldmvnorm(obs = obs[p,], chol = aperm(C, perm = p))
ll3 <- ldmvnorm(obs = obs, invchol = L)
ll4 <- ldmvnorm(obs = obs[p,], invchol = aperm(L, perm = p))
chk(ll1, ll2)
chk(ll1, ll3)
chk(ll1, ll4)

### C
### diag = TRUE w/o stand
ll <- function(x) {
    C <- as.chol(ltMatrices(x, diag = TRUE, byrow = BYROW))
    Ct <- aperm(C, perm = p)
    -ldmvnorm(obs = obs[p,], chol = Ct)
}

s <- function(x) {
    C <- as.chol(ltMatrices(x, diag = TRUE, byrow = BYROW))
    Ct <- aperm(C, perm = p)
    sC <- sldmvnorm(obs = obs[p,], chol = Ct)$chol
    ret <- deperma(chol = C, permuted_chol = Ct, perm = p, score_schol = sC)
    -rowSums(Lower_tri(ret, diag = TRUE))
}

g1 <- grad(ll, c(C))
s1 <- s(c(C))
chk(g1, s1)

op1 <- optim(c(C), fn = ll, gr = s, method = "L-BFGS-B")
max(abs(ltMatrices(op1$par, diag = TRUE, byrow = BYROW) - C))

### check against unpermuted (expect same results)
ll <- function(x) {
    C <- ltMatrices(x, diag = TRUE, byrow = BYROW)
    -ldmvnorm(obs = obs, chol = C)
}

s <- function(x) {
    C <- ltMatrices(x, diag = TRUE, byrow = BYROW)
    ret <- sldmvnorm(obs = obs, chol = C)$chol
    -rowSums(Lower_tri(ret, diag = TRUE))
}

op2 <- optim(c(C), fn = ll, gr = s, method = "L-BFGS-B")
chk(max(abs(ltMatrices(op2$par, diag = TRUE, byrow = BYROW) - C)) < tol, TRUE)

chk(op1, op2)

### diag = FALSE
Cd <- ltMatrices(runif(J * (J - 1) / 2), byrow = BYROW)

### w/ standardisation (1. stand, 2. perm)
ll <- function(x) {
    C <- as.chol(ltMatrices(x, diag = FALSE, byrow = BYROW))
    Cs <- standardize(chol = C)
    Ct <- aperm(Cs, perm = p)
    -ldmvnorm(obs = obs[p,], chol = Ct)
}

s <- function(x) {
    C <- ltMatrices(x, diag = FALSE, byrow = BYROW)
    Cs <- standardize(chol = C)
    Ct <- aperm(Cs, perm = p)
    sC <- sldmvnorm(obs = obs[p,], chol = Ct)$chol
    ret <- deperma(chol = Cs, permuted_chol = Ct, perm = p, score_schol = sC)
    ret <- destandardize(chol = C, score_schol = ret)
    -rowSums(Lower_tri(ret, diag = FALSE))
}

chk(grad(ll, c(Cd)), s(c(Cd)))

### w/o standardisation 
ll <- function(x) {
    C <- as.chol(ltMatrices(x, diag = FALSE, byrow = BYROW))
    Ct <- aperm(C, perm = p)
    -ldmvnorm(obs = obs[p,], chol = Ct)
}

s <- function(x) {
    C <- as.chol(ltMatrices(x, diag = FALSE, byrow = BYROW))
    diagonals(C) <- 1		### deperma expects diagonals
    Ct <- aperm(as.chol(C), perm = p)
    sC <- sldmvnorm(obs = obs[p,], chol = Ct)$chol
    ret <- deperma(chol = C, permuted_chol = Ct, perm = p, score_schol = sC)
    -rowSums(Lower_tri(ret, diag = FALSE))
}

chk(grad(ll, c(Cd)), s(c(Cd)))

### L
### diag = TRUE w/o stand
ll <- function(x) {
    C <- as.invchol(ltMatrices(x, diag = TRUE, byrow = BYROW))
    Ct <- aperm(C, perm = p)
    -ldmvnorm(obs = obs[p,], invchol = Ct)
}

s <- function(x) {
    C <- as.invchol(ltMatrices(x, diag = TRUE, byrow = BYROW))
    Ct <- aperm(C, perm = p)
    sC <- sldmvnorm(obs = obs[p,], invchol = Ct)$invchol
    ret <- deperma(invchol = C, permuted_invchol = Ct, perm = p, score_schol = -vectrick(Ct, sC))
    -rowSums(Lower_tri(ret, diag = TRUE))
}

g2 <- grad(ll, c(L))
chk(g2, s(c(L)))
chk(g2, c(Lower_tri(-vectrick(C, ltMatrices(g1, byrow = BYROW, diag = TRUE)), diag = TRUE)))

op3 <- optim(c(L), fn = ll, gr = s, method = "L-BFGS-B")
chk(max(abs(ltMatrices(op3$par, diag = TRUE, byrow = BYROW) - L)) < tol, TRUE)

### check against unpermuted (expect same results)
ll <- function(x) {
    C <- ltMatrices(x, diag = TRUE, byrow = BYROW)
    -ldmvnorm(obs = obs, invchol = C)
}

s <- function(x) {
    C <- ltMatrices(x, diag = TRUE, byrow = BYROW)
    ret <- sldmvnorm(obs = obs, invchol = C)$invchol
    -rowSums(Lower_tri(ret, diag = TRUE))
}

op4 <- optim(c(L), fn = ll, gr = s, method = "L-BFGS-B")
chk(max(abs(ltMatrices(op4$par, diag = TRUE, byrow = BYROW) - L)) < tol, TRUE)

### diag = FALSE
Ld <- ltMatrices(runif(J * (J - 1) / 2), byrow = BYROW)

### w/ standardisation (1. stand, 2. perm)
ll <- function(x) {
    C <- as.invchol(ltMatrices(x, diag = FALSE, byrow = BYROW))
    Cs <- standardize(invchol = C)
    Ct <- aperm(Cs, perm = p)
    -ldmvnorm(obs = obs[p,], invchol = Ct)
}

s <- function(x) {
    C <- as.invchol(ltMatrices(x, diag = FALSE, byrow = BYROW))
    Cs <- standardize(invchol = C)
    Ct <- aperm(Cs, perm = p)
    sC <- sldmvnorm(obs = obs[p,], invchol = Ct)$invchol
    ret <- deperma(invchol = Cs, permuted_invchol = Ct, perm = p, score_schol = -vectrick(Ct, sC))
    ret <- destandardize(invchol = C, score_schol = -vectrick(Cs, ret))
    -rowSums(Lower_tri(ret, diag = FALSE))
}

chk(grad(ll, c(Ld)), s(c(Ld)))

### w/o standardisation 
ll <- function(x) {
    C <- as.invchol(ltMatrices(x, diag = FALSE, byrow = BYROW))
    Ct <- aperm(C, perm = p)
    -ldmvnorm(obs = obs[p,], invchol = Ct)
}

s <- function(x) {
    C <- as.invchol(ltMatrices(x, diag = FALSE, byrow = BYROW))
    diagonals(C) <- 1		### deperma expects diagonals
    Ct <- aperm(as.invchol(C), perm = p)
    sC <- sldmvnorm(obs = obs[p,], invchol = Ct)$invchol
    ret <- deperma(invchol = C, permuted_invchol = Ct, perm = p, score_schol = -vectrick(Ct, sC))
    -rowSums(Lower_tri(ret, diag = FALSE))
}

chk(grad(ll, c(Ld)), s(c(Ld)))
})

BYROW <- FALSE
EVAL(thischeck)

BYROW <- TRUE
EVAL(thischeck)
