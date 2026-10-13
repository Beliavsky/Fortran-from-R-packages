# SPDX-License-Identifier: GPL-2.0-or-later
# Run from the package root:
# Rscript validation/generate_design_oracle.R /path/to/CRAN/irr [output.txt]
# Requires the original lpSolve R package for N2.cohen.kappa.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("Supply the unpacked upstream irr directory")
upstream <- args[1]
if (length(args) >= 2L) sink(args[2])
options(digits = 17, width = 110)
cat("R oracle:", R.version.string, "\n")
for (file in c("N.cohen.kappa.R", "N2.cohen.kappa.R", "relInterIntra.R")) {
  source(file.path(upstream, "R", file))
}
cat("Upstream irr version:", read.dcf(file.path(upstream, "DESCRIPTION"))[1,"Version"], "\n")
cat("N.cohen.kappa\n")
dput(c(N.cohen.kappa(.5, .5, .85, .7),
       N.cohen.kappa(.3, .4, .55, .25, twosided = TRUE),
       N.cohen.kappa(.5, .5, -.2, -.5, alpha = .01, power = .9),
       N.cohen.kappa(.5, .5, 1, .7)))
if (requireNamespace("lpSolve", quietly = TRUE)) {
  lp <- lpSolve::lp
  cat("N2.cohen.kappa, lpSolve", as.character(packageVersion("lpSolve")), "\n")
  dput(c(N2.cohen.kappa(c(.2, .25, .55), .6, .4),
         N2.cohen.kappa(c(.2, .25, .55), .6, .4, power = .9),
         N2.cohen.kappa(c(.2, .05, .2, .05, .2, .3), .5, .1),
         N2.cohen.kappa(c(.5, .5), .85, .7),
         N2.cohen.kappa(rep(.25, 4), .7, .4, twosided = TRUE),
         N2.cohen.kappa(rep(.1, 10), .6, .3),
         N2.cohen.kappa(c(0, .3, .7), .7, .4)))
} else {
  cat("N2.cohen.kappa skipped: install upstream lpSolve to regenerate all oracles.\n")
}
table4 <- matrix(c(
 -2,16,5,11,7,-7,18,4,0,0,-3,3,7,-6,1,-13,2,4,-10,8,7,-3,-5,5,0,7,-8,1,-3,
 0,16,6,10,8,-8,19,5,-3,0,-2,-1,9,-7,1,-14,1,4,-9,9,6,-2,-5,5,-1,6,-8,1,-3,
 1,15,6,10,6,-8,19,5,-2,-2,-2,1,9,-6,0,-14,0,3,-10,8,7,-4,-7,5,-1,6,-8,2,-3,
 2,12,4,9,5,-9,17,5,-7,1,-4,-1,4,-8,-2,-12,-1,7,-10,2,8,-5,-6,3,-4,4,-10,1,-5,
 1,14,4,7,6,-10,17,5,-6,2,-3,-2,4,-10,-2,-12,0,6,-11,8,7,-5,-8,4,-3,4,-11,-1,-4,
 1,13,4,8,6,-9,17,5,-5,1,-3,1,2,-9,-3,-12,0,4,-10,8,7,-5,-7,4,-4,4,-10,0,-5), ncol = 6)
cat("relInterIntra goniometer example, all five numeric vectors\n")
dput(relInterIntra(table4, nrater = 2)$value)
x <- outer(1:7, 1:6, function(i, j) .4*i^2 + .3*j + ((i*j+j^2) %% 7)*.7)
cat("relInterIntra deterministic 7 subjects, 3 raters, 2 repetitions, custom options\n")
dput(relInterIntra(x, nrater = 3, rho0inter = .2, rho0intra = .3, conf.level = .9)$value)
if (length(args) >= 2L) sink()
