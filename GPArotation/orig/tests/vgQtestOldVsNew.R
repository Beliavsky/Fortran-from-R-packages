# vgQtestOldVsNew.R
# Regression tests: new vgQ functions vs original implementations.
# f and Gq are compared within tolerance to allow for floating point
# differences due to different order of arithmetic operations.
# Method strings are not compared --- changes are intentional improvements.
# A real bug (mccammon Gq) was detected by these tests in 2026.4-1.

require("GPArotation")
source("vgQOLD.R")

all.ok <- TRUE
tol    <- 1e-10

set.seed(123)
L8x2  <- matrix(rnorm(16),  nrow = 8,  ncol = 2)
L10x3 <- matrix(rnorm(30),  nrow = 10, ncol = 3)
L15x4 <- matrix(rnorm(60),  nrow = 15, ncol = 4)
L20x5 <- matrix(rnorm(100), nrow = 20, ncol = 5)

matrices <- list(L8x2, L10x3, L15x4, L20x5)

# Helper: check structural contract
check_contract <- function(res, name) {
  ok <- TRUE
  if (!all(c("f", "Gq", "Method") %in% names(res))) {
    cat("FAILED:", name, "- missing list elements\n")
    ok <- FALSE
  }
  if (!is.numeric(res$f) || length(res$f) != 1) {
    cat("FAILED:", name, "- f is not a scalar numeric\n")
    ok <- FALSE
  }
  if (!is.character(res$Method) || length(res$Method) != 1) {
    cat("FAILED:", name, "- Method is not a scalar character\n")
    ok <- FALSE
  }
  ok
}

# Helper: check new vs old within tolerance
# Method strings are not compared --- changes are intentional
check_values <- function(new, old, name, tol = 1e-10) {
  ok <- TRUE
  if (max(abs(new$f - old$f)) > tol) {
    cat("FAILED:", name, "- f differs beyond tolerance\n")
    cat("  max diff:", max(abs(new$f - old$f)), "\n")
    ok <- FALSE
  }
  if (max(abs(new$Gq - old$Gq)) > tol) {
    cat("FAILED:", name, "- Gq differs beyond tolerance\n")
    cat("  max diff:", max(abs(new$Gq - old$Gq)), "\n")
    ok <- FALSE
  }
  ok
}

# -------------------------------------------------------------------
# oblimin --- gam = 0, 0.25, 0.5
# -------------------------------------------------------------------
for (L in matrices) {
  for (gam in c(0, 0.25, 0.5)) {
    nm  <- paste("vgQ.oblimin gam=", gam, dim(L)[1], "x", dim(L)[2])
    new <- GPArotation:::vgQ.oblimin(L, gam = gam)
    old <- vgQOLD.oblimin(L, gam = gam)
    all.ok <- all.ok & check_contract(new, nm)
    all.ok <- all.ok & check_values(new, old, nm)
  }
}

# -------------------------------------------------------------------
# quartimin
# -------------------------------------------------------------------
for (L in matrices) {
  nm  <- paste("vgQ.quartimin", dim(L)[1], "x", dim(L)[2])
  new <- GPArotation:::vgQ.quartimin(L)
  old <- vgQOLD.quartimin(L)
  all.ok <- all.ok & check_contract(new, nm)
  all.ok <- all.ok & check_values(new, old, nm)
}

# -------------------------------------------------------------------
# target --- with and without NA
# -------------------------------------------------------------------
set.seed(456)
for (L in matrices) {
  Target <- matrix(rnorm(length(L)), nrow = nrow(L), ncol = ncol(L))
  nm  <- paste("vgQ.target", dim(L)[1], "x", dim(L)[2])
  new <- GPArotation:::vgQ.target(L, Target = Target)
  old <- vgQOLD.target(L, Target = Target)
  all.ok <- all.ok & check_contract(new, nm)
  all.ok <- all.ok & check_values(new, old, nm)
}

# target with NA
Target_na <- matrix(c(rep(NA, 4), rep(0, 4)), nrow = 8, ncol = 2)
new <- GPArotation:::vgQ.target(L8x2, Target = Target_na)
old <- vgQOLD.target(L8x2, Target = Target_na)
all.ok <- all.ok & check_contract(new, "vgQ.target with NA")
all.ok <- all.ok & check_values(new, old, "vgQ.target with NA")

# -------------------------------------------------------------------
# pst
# -------------------------------------------------------------------
set.seed(789)
for (L in matrices) {
  W      <- matrix(sample(0:1, length(L), replace = TRUE),
                   nrow = nrow(L), ncol = ncol(L))
  Target <- matrix(0, nrow = nrow(L), ncol = ncol(L))
  nm  <- paste("vgQ.pst", dim(L)[1], "x", dim(L)[2])
  new <- GPArotation:::vgQ.pst(L, W = W, Target = Target)
  old <- vgQOLD.pst(L, W = W, Target = Target)
  all.ok <- all.ok & check_contract(new, nm)
  all.ok <- all.ok & check_values(new, old, nm)
}

# -------------------------------------------------------------------
# oblimax
# -------------------------------------------------------------------
for (L in matrices) {
  nm  <- paste("vgQ.oblimax", dim(L)[1], "x", dim(L)[2])
  new <- GPArotation:::vgQ.oblimax(L)
  old <- vgQOLD.oblimax(L)
  all.ok <- all.ok & check_contract(new, nm)
  all.ok <- all.ok & check_values(new, old, nm)
}

# -------------------------------------------------------------------
# entropy
# -------------------------------------------------------------------
for (L in matrices) {
  nm  <- paste("vgQ.entropy", dim(L)[1], "x", dim(L)[2])
  new <- GPArotation:::vgQ.entropy(L)
  old <- vgQOLD.entropy(L)
  all.ok <- all.ok & check_contract(new, nm)
  all.ok <- all.ok & check_values(new, old, nm)
}

# -------------------------------------------------------------------
# quartimax
# -------------------------------------------------------------------
for (L in matrices) {
  nm  <- paste("vgQ.quartimax", dim(L)[1], "x", dim(L)[2])
  new <- GPArotation:::vgQ.quartimax(L)
  old <- vgQOLD.quartimax(L)
  all.ok <- all.ok & check_contract(new, nm)
  all.ok <- all.ok & check_values(new, old, nm)
}

# -------------------------------------------------------------------
# varimax
# -------------------------------------------------------------------
for (L in matrices) {
  nm  <- paste("vgQ.varimax", dim(L)[1], "x", dim(L)[2])
  new <- GPArotation:::vgQ.varimax(L)
  old <- vgQOLD.varimax(L)
  all.ok <- all.ok & check_contract(new, nm)
  all.ok <- all.ok & check_values(new, old, nm)
}

# -------------------------------------------------------------------
# simplimax --- default k and custom k
# -------------------------------------------------------------------
for (L in matrices) {
  for (k in c(nrow(L), floor(nrow(L) / 2))) {
    nm  <- paste("vgQ.simplimax k=", k, dim(L)[1], "x", dim(L)[2])
    new <- GPArotation:::vgQ.simplimax(L, k = k)
    old <- vgQOLD.simplimax(L, k = k)
    all.ok <- all.ok & check_contract(new, nm)
    all.ok <- all.ok & check_values(new, old, nm)
  }
}

# -------------------------------------------------------------------
# bentler
# -------------------------------------------------------------------
for (L in matrices) {
  nm  <- paste("vgQ.bentler", dim(L)[1], "x", dim(L)[2])
  new <- GPArotation:::vgQ.bentler(L)
  old <- vgQOLD.bentler(L)
  all.ok <- all.ok & check_contract(new, nm)
  all.ok <- all.ok & check_values(new, old, nm)
}

# -------------------------------------------------------------------
# tandemI
# -------------------------------------------------------------------
for (L in matrices) {
  nm  <- paste("vgQ.tandemI", dim(L)[1], "x", dim(L)[2])
  new <- GPArotation:::vgQ.tandemI(L)
  old <- vgQOLD.tandemI(L)
  all.ok <- all.ok & check_contract(new, nm)
  all.ok <- all.ok & check_values(new, old, nm)
}

# -------------------------------------------------------------------
# tandemII
# -------------------------------------------------------------------
for (L in matrices) {
  nm  <- paste("vgQ.tandemII", dim(L)[1], "x", dim(L)[2])
  new <- GPArotation:::vgQ.tandemII(L)
  old <- vgQOLD.tandemII(L)
  all.ok <- all.ok & check_contract(new, nm)
  all.ok <- all.ok & check_values(new, old, nm)
}

# -------------------------------------------------------------------
# geomin --- default and custom delta
# -------------------------------------------------------------------
for (L in matrices) {
  for (delta in c(0.01, 0.001)) {
    nm  <- paste("vgQ.geomin delta=", delta, dim(L)[1], "x", dim(L)[2])
    new <- GPArotation:::vgQ.geomin(L, delta = delta)
    old <- vgQOLD.geomin(L, delta = delta)
    all.ok <- all.ok & check_contract(new, nm)
    all.ok <- all.ok & check_values(new, old, nm)
  }
}

# -------------------------------------------------------------------
# bigeomin
# -------------------------------------------------------------------
for (L in matrices) {
  nm  <- paste("vgQ.bigeomin", dim(L)[1], "x", dim(L)[2])
  new <- GPArotation:::vgQ.bigeomin(L)
  old <- vgQOLD.bigeomin(L)
  all.ok <- all.ok & check_contract(new, nm)
  all.ok <- all.ok & check_values(new, old, nm)
}

# -------------------------------------------------------------------
# cf --- kappa = 0, 0.3, 1
# -------------------------------------------------------------------
for (L in matrices) {
  for (kappa in c(0, 0.3, 1)) {
    nm  <- paste("vgQ.cf kappa=", kappa, dim(L)[1], "x", dim(L)[2])
    new <- GPArotation:::vgQ.cf(L, kappa = kappa)
    old <- vgQOLD.cf(L, kappa = kappa)
    all.ok <- all.ok & check_contract(new, nm)
    all.ok <- all.ok & check_values(new, old, nm)
  }
}

# -------------------------------------------------------------------
# infomax
# -------------------------------------------------------------------
for (L in matrices) {
  nm  <- paste("vgQ.infomax", dim(L)[1], "x", dim(L)[2])
  new <- GPArotation:::vgQ.infomax(L)
  old <- vgQOLD.infomax(L)
  all.ok <- all.ok & check_contract(new, nm)
  all.ok <- all.ok & check_values(new, old, nm)
}

# -------------------------------------------------------------------
# mccammon -- this test no longer valid as of May 6, 2026
# -------------------------------------------------------------------
#for (L in matrices) {
#  nm  <- paste("vgQ.mccammon", dim(L)[1], "x", dim(L)[2])
#  new <- GPArotation:::vgQ.mccammon(L)
#  old <- vgQOLD.mccammon(L)
#  all.ok <- all.ok & check_contract(new, nm)
#  all.ok <- all.ok & check_values(new, old, nm)
#}

# -------------------------------------------------------------------
# bifactor --- first column of Gq must be zero
# -------------------------------------------------------------------
for (L in matrices) {
  nm  <- paste("vgQ.bifactor", dim(L)[1], "x", dim(L)[2])
  new <- GPArotation:::vgQ.bifactor(L)
  old <- vgQOLD.bifactor(L)
  all.ok <- all.ok & check_contract(new, nm)
  all.ok <- all.ok & check_values(new, old, nm)
  if (!all(new$Gq[, 1] == 0)) {
    cat("FAILED:", nm, "- first column of Gq not zero\n")
    all.ok <- FALSE
  }
}

# -------------------------------------------------------------------
# varimin
# -------------------------------------------------------------------
for (L in matrices) {
  nm  <- paste("vgQ.varimin", dim(L)[1], "x", dim(L)[2])
  new <- GPArotation:::vgQ.varimin(L)
  old <- vgQOLD.varimin(L)
  all.ok <- all.ok & check_contract(new, nm)
  all.ok <- all.ok & check_values(new, old, nm)
}

# -------------------------------------------------------------------
# lp.wls
# -------------------------------------------------------------------
set.seed(321)
for (L in matrices) {
  W   <- matrix(runif(length(L)), nrow = nrow(L), ncol = ncol(L))
  nm  <- paste("vgQ.lp.wls", dim(L)[1], "x", dim(L)[2])
  new <- GPArotation:::vgQ.lp.wls(L, W = W)
  old <- vgQOLD.lp.wls(L, W = W)
  all.ok <- all.ok & check_contract(new, nm)
  all.ok <- all.ok & check_values(new, old, nm)
}

# -------------------------------------------------------------------
# Gq dimensions match L for all criteria
# -------------------------------------------------------------------
criteria <- c("oblimin", "quartimin", "oblimax", "entropy", "quartimax",
              "varimax", "simplimax", "bentler", "tandemI", "tandemII",
              "geomin", "bigeomin", "infomax", "mccammon", "bifactor",
              "varimin")

for (crit in criteria) {
  fn <- get(paste("vgQ", crit, sep = "."),
            envir = asNamespace("GPArotation"))
  for (L in matrices) {
    res <- fn(L)
    if (!identical(dim(res$Gq), dim(L))) {
      cat("FAILED: vgQ.", crit, "Gq dimensions do not match L:",
          dim(L), "vs", dim(res$Gq), "\n")
      all.ok <- FALSE
    }
  }
}

cat("vgQ regression tests completed.\n")
if (!all.ok) stop("some tests FAILED")