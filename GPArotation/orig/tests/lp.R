require("GPArotation")

# Local sort function --- used for lp tests since lpT/lpQ are not
# standard GPA rotations and require local sorting for consistency.
sortFac <- function(x) {
  vx   <- order(colSums(x$loadings^2), decreasing = TRUE)
  Dsgn <- diag(sign(colSums(x$loadings^3)))[, vx]
  x$Th      <- x$Th %*% Dsgn
  x$loadings <- x$loadings %*% Dsgn
  if ("Phi" %in% names(x))
    x$Phi <- diag(1 / diag(Dsgn)) %*% x$Phi %*% Dsgn
  x
}

fuzz   <- 1e-5  # using eps=1e-5 tests do not do better than this
all.ok <- TRUE


# --- Synthetic 30x3 data ---
L <- rbind(diag(3), diag(3), diag(3), diag(3), diag(3),
           diag(3), diag(3), diag(3), diag(3), diag(3))
True_rot <- matrix(c(1, 0.02079577, 0.5024378,
                     0, 0.99978374, 0.2635086,
                     0, 0,          0.8234801),
                   3, 3, byrow = TRUE)
L1 <- L %*% t(True_rot)

# Repeating pattern helper
rep3 <- function(r1, r2, r3, n = 10)
  t(matrix(rep(c(r1, r2, r3), n), 3, 3 * n))


# --- Test 1: lpQ default p, fixed identity start ---

r1 <- lpQ(L1, Tmat = diag(3), maxit = 1000)

tst1 <- rep3(
  c(6.077051499e-05, 9.990889331e-01, 0.0018179274126),
  c(9.997552961e-01, 6.087759033e-05, 0.0008923378188),
  c(8.931676374e-04, 1.821872386e-03, 0.9988446821739))

if (fuzz < max(abs(sortFac(r1)$loadings - tst1))) {
  cat("Test 1 lpQ default p: value not same as test value.\n")
  print(unclass(r1$loadings), digits = 18)
  cat("difference:\n")
  print(unclass(sortFac(r1)$loadings) - tst1, digits = 18)
  all.ok <- FALSE
}


# --- Test 2: lpQ p=0.5, fixed start ---
# p=0.5 gives sparser solution than default --- genuinely different from test 1

r2 <- lpQ(L1, p = 0.5, Tmat = diag(3), maxit = 1000)

tst2 <- rep3(
  c(3.681422987e-06, 9.999550303e-01, 8.936779409e-05),
  c(9.999865823e-01, 3.696883647e-06, 4.870092977e-05),
  c(4.870461564e-05, 8.935912376e-05, 9.999417744e-01))

if (fuzz < max(abs(sortFac(r2)$loadings - tst2))) {
  cat("Test 2 lpQ p=0.5: value not same as test value.\n")
  print(unclass(r2$loadings), digits = 18)
  cat("difference:\n")
  print(unclass(sortFac(r2)$loadings) - tst2, digits = 18)
  all.ok <- FALSE
}


# --- Test 3: lpQ on NetherlandsTV 2-factor, fixed start ---

data("WansbeekMeijer", package = "GPArotation")
fa.nl2 <- factanal(factors = 2, covmat = NetherlandsTV,
                   rotation = "none")

r3 <- lpQ(fa.nl2$loadings, p = 0.75,
          Tmat = diag(2), maxit = 1000)

tst3 <- matrix(c(
  -0.002173520829,  0.792259805380,
   0.100607661633,  0.781466305964,
   0.002154155542,  0.772058958326,
   0.617378882323,  0.138659737124,
   0.707273557353,  0.090568574748,
   0.822844605470, -0.009239665513,
   0.725226973378,  0.002693412417
), ncol = 2, nrow = 7, byrow = TRUE)

if (fuzz < max(abs(sortFac(r3)$loadings - tst3))) {
  cat("Test 3 lpQ NetherlandsTV: value not same as test value.\n")
  print(unclass(r3$loadings), digits = 18)
  cat("difference:\n")
  print(unclass(sortFac(r3)$loadings) - tst3, digits = 18)
  all.ok <- FALSE
}


# --- Test 4: lpT on NetherlandsTV 3-factor, fixed start ---

fa.nl3 <- factanal(factors = 3, covmat = NetherlandsTV,
                   rotation = "none")

r4 <- lpT(fa.nl3$loadings, p = 0.75,
          Tmat = diag(3), maxit = 1000)

tst4 <- matrix(c(
   0.3175272930,  0.7249593608,  -0.0190811964,
   0.3837479919,  0.7455368516,   0.0506357055,
   0.3015091348,  0.7115867577,   0.0015401560,
   0.4703884148,  0.3165308093,   0.4592056754,
   0.5706091490,  0.2747953193,   0.4004865243,
   0.5997665911,  0.2228087629,   0.4986242377,
   0.9974965731,  0.0003420514,  -0.0006889090
), ncol = 3, nrow = 7, byrow = TRUE)

if (fuzz < max(abs(sortFac(r4)$loadings - tst4))) {
  cat("Test 4 lpT NetherlandsTV: value not same as test value.\n")
  print(unclass(r4$loadings), digits = 18)
  cat("difference:\n")
  print(unclass(sortFac(r4)$loadings) - tst4, digits = 18)
  all.ok <- FALSE
}


if (!all.ok) stop("some tests FAILED")