# Reference calculations for compositions-fortran v0.1.0
# Run with R after installing the source package compositions.
library(compositions)

x <- acomp(c(2, 3, 5))
cat("clr:\n"); print(unclass(clr(x)))
cat("ilr:\n"); print(unclass(ilr(x)))
cat("alr:\n"); print(unclass(alr(x)))
cat("apt:\n"); print(unclass(apt(rcomp(x))))

set.seed(12345)
xd <- rDirichlet.acomp(10000, c(2,3,5))
cat("Dirichlet sample mean:\n"); print(colMeans(unclass(xd)))
cat("Dirichlet fit:\n"); print(fitDirichlet(xd)$alpha)

z <- acomp(rbind(c(.6,.3,.1), c(.5,.35,.15), c(.4,.4,.2)))
cat("variation:\n"); print(variation(z))
cat("PBmaxvar basis:\n"); print(gsi.PrinBal(z, method="PBmaxvar"))
