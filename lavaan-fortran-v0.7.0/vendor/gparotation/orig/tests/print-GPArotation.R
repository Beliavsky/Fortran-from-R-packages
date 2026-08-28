# Tests for random starts and factor ordering/sorting.
# Verifies that:
#   1. Different random starts produce different raw output (before sorting)
#   2. Sorting produces identical output regardless of starting point
#   3. factanal built-in rotation agrees with GPArotation two-step (R >= 4.5.1)

require("GPArotation")

require("stats")
require("GPArotation")

fuzz   <- 1e-5
all.ok <- TRUE

athl <- matrix(c(
   .73, -.07,  .50,  .82, -.01,  .27,  .77, -.46, -.22,
   .78,  .17,  .03,  .77,  .41,  .13,  .81, -.01,  .27,
   .71, -.45, -.30,  .82,  .12, -.11,  .66, -.15, -.45,
   .39,  .76, -.40),
  byrow = TRUE, ncol = 3)

# Seeds 238 and 46 were chosen because they converge to the same local
# minimum, allowing us to test that sorting produces identical output
# regardless of starting point.
set.seed(238)
z1 <- quartimin(athl, Tmat = Random.Start(3), algorithm = "legacy", fwindow = 1)

set.seed(46)
z2 <- quartimin(athl, Tmat = Random.Start(3), algorithm = "legacy", fwindow = 1)

# --- Before sorting: z1 and z2 should differ ---
if (fuzz > max(abs(z1$loadings - z2$loadings))) {
  cat("Random starts test failed: loadings should differ before sorting\n")
  all.ok <- FALSE
}

if (fuzz > max(abs(z1$Th - z2$Th))) {
  cat("Random starts test failed: Th should differ before sorting\n")
  all.ok <- FALSE
}

if (fuzz > max(abs(z1$Phi - z2$Phi))) {
  cat("Random starts test failed: Phi should differ before sorting\n")
  all.ok <- FALSE
}

if (!all.ok) stop("some tests FAILED before sorting")

# --- After sorting: z1 and z2 should agree ---
z1s <- print(z1, sortLoadings = TRUE, rotateMat = TRUE, Table = TRUE)
z2s <- print(z2, sortLoadings = TRUE, rotateMat = TRUE, Table = TRUE)

all.ok <- TRUE

if (fuzz < max(abs(z1s$loadings - z2s$loadings))) {
  cat("Sorting test failed: loadings should agree after sorting\n")
  cat("difference:\n")
  print(z1s$loadings - z2s$loadings, digits = 18)
  all.ok <- FALSE
}

if (fuzz < max(abs(z1s$Th - z2s$Th))) {
  cat("Sorting test failed: Th should agree after sorting\n")
  cat("difference:\n")
  print(z1s$Th - z2s$Th, digits = 18)
  all.ok <- FALSE
}

if (fuzz < max(abs(z1s$Phi - z2s$Phi))) {
  cat("Sorting test failed: Phi should agree after sorting\n")
  cat("difference:\n")
  print(z1s$Phi - z2s$Phi, digits = 18)
  all.ok <- FALSE
}

if (!all.ok) stop("some tests FAILED after sorting")

# --- factanal regression test (R >= 4.5.1) ---
# Prior to R 4.5.1, factanal had a bug in factor reordering after rotation.
# This was reported by Bernaards and others and fixed by the R core team.
# The two-step procedure (factanal rotation="none" + GPArotation) is always
# correct regardless of R version.

data(ability.cov)
set.seed(134)
Tmat <- Random.Start(3)

fa_unrotated <- factanal(factors = 3, covmat = ability.cov, rotation = "none")
gpa          <- oblimin(loadings(fa_unrotated), Tmat = Tmat)
gpa_sorted   <- print(gpa, sortLoadings = TRUE)

if (getRversion() >= "4.5.1") {
  set.seed(134)
  fa_rotated <- factanal(factors = 3, covmat = ability.cov, rotation = "oblimin",
                         control = list(rotate = list(Tmat = Tmat)))

  # Use abs() to account for differing sign correction heuristics between
  # base R and GPArotation
  if (1e-4 < max(abs(abs(fa_rotated$loadings) - abs(gpa_sorted$loadings)))) {
    cat("factanal regression test failed: loadings disagree\n")
    cat("difference:\n")
    print(abs(fa_rotated$loadings) - abs(gpa_sorted$loadings), digits = 18)
    all.ok <- FALSE
  }

  if (1e-4 < max(abs(fa_rotated$Phi - gpa_sorted$Phi))) {
    cat("factanal regression test failed: Phi disagrees\n")
    cat("difference:\n")
    print(fa_rotated$Phi - gpa_sorted$Phi, digits = 18)
    all.ok <- FALSE
  }

} else {
  cat("Skipping factanal ordering test — requires R >= 4.5.1\n")
  cat("Use the two-step procedure (factanal rotation='none' + GPArotation) for correct results.\n")
}

cat("tests completed.\n")

if (!all.ok) stop("some tests FAILED")

#> z1
#Oblique rotation method Quartimin converged.
#Loadings:
#         [,1]    [,2]     [,3]
# [1,]  0.9451 -0.0535 -0.18033
# [2,]  0.7725  0.1431  0.01187
# [3,]  0.1323  0.8617 -0.12847
# [4,]  0.5377  0.2050  0.28753
# [5,]  0.6888 -0.0555  0.44072
# [6,]  0.7665  0.1386  0.00987
# [7,]  0.0150  0.8967 -0.08931
# [8,]  0.4047  0.3792  0.32647
# [9,] -0.1056  0.7915  0.24071
#[10,] -0.0155 -0.0165  0.94994
#
#                [,1]  [,2]  [,3]
#SS loadings    3.034 2.405 1.401
#Proportion Var 0.303 0.240 0.140
#Cumulative Var 0.303 0.544 0.684
#
#Phi:
#      [,1]  [,2]  [,3]
#[1,] 1.000 0.554 0.259
#[2,] 0.554 1.000 0.186
#[3,] 0.259 0.186 1.000
#> z2
#Oblique rotation method Quartimin converged.
#Loadings:
#         [,1]    [,2]     [,3]
# [1,]  0.9451 -0.0535 -0.18033
# [2,]  0.7725  0.1431  0.01187
# [3,]  0.1323  0.8617 -0.12847
# [4,]  0.5377  0.2050  0.28753
# [5,]  0.6888 -0.0555  0.44072
# [6,]  0.7665  0.1386  0.00987
# [7,]  0.0150  0.8967 -0.08930
# [8,]  0.4047  0.3792  0.32647
# [9,] -0.1056  0.7915  0.24071
#[10,] -0.0155 -0.0165  0.94994
#
#                [,1]  [,2]  [,3]
#SS loadings    3.034 2.405 1.401
#Proportion Var 0.303 0.240 0.140
#Cumulative Var 0.303 0.544 0.684
#
#Phi:
#      [,1]  [,2]  [,3]
#[1,] 1.000 0.554 0.259
#[2,] 0.554 1.000 0.186
#[3,] 0.259 0.186 1.000


##########################################################
# RUNNING A PRINT WITHOUT ERRORS
##########################################################
#  data(ability.cov)
#  L <- loadings(factanal(factors = 2, covmat=ability.cov))
#
#v <- GPFRSoblq(L, eps = 1e-7, method = "oblimin", methodArgs = list(gam = .5), randomStarts = 100)
#GPArotation:::print.GPArotation(v, rotateMat = T, Table = T)
#print(v, rotateMat = T, Table = T)
#
#GPArotation:::summary.GPArotation(v) 
#summary(v)


