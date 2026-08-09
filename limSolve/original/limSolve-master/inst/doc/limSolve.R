## ----setup, include=FALSE, warning=FALSE, message=FALSE-----------------------
knitr::opts_chunk$set(echo = TRUE)
require(limSolve)

## -----------------------------------------------------------------------------
E <- matrix(nrow = 3, ncol = 3,
            data = c(3, 1, 2, 2, 0, 0, 1, 0, 2))
F <- c(2, 1, 8)

solve(E, F)
Solve(E, F)

## -----------------------------------------------------------------------------
E2 <- rbind(E, E[1, ] + E[2, ])
F2 <- c(F, F[1] + F[2])

#solve(E2,F2) # error 'a' (4 x 3) must be square
Solve(E2, F2)

## -----------------------------------------------------------------------------
E <- matrix(nrow = 4, ncol = 3,
            data = c(3, 1, 2, 0, 2, 0, 0, 1, 1, 0, 2, 0))
F <- c(2, 1, 8, 3)

lsei(A = E, B = F, fulloutput = TRUE, verbose = FALSE)

## -----------------------------------------------------------------------------
Solve(E, F)

## -----------------------------------------------------------------------------
G <- matrix(nrow = 2, ncol = 3, byrow = TRUE,
            data = c(-1, 2, 0, 1, 0, -1))
H <- c(-3, 2)

lsei(E = E, F = F, G = G, H = H, verbose = FALSE)

## -----------------------------------------------------------------------------
lsei(A = E, B = F, G = G, H = H, verbose = FALSE)

## -----------------------------------------------------------------------------
E <- matrix(nrow = 2, ncol = 3,
            data = c(3, 1, 2, 0, 1, 0))
F <- c(2, 1)

Solve(E, F)
ldei(E = E, F = F)$X

## -----------------------------------------------------------------------------
lsei(E = E, F = F, A = diag(3), B = rep(0, 3), verbose = FALSE)$X

## ----u1b, fig.cap="Fig: Random sample of the underdetermined system including only equalities."----
xs <- xsample(E = E, F = F, iter = 500)$X
plot(xs[ ,2], xs[ ,3])

## -----------------------------------------------------------------------------
E <- matrix(ncol = 4, byrow = TRUE,
            data = c(3, 2, 1, 4, 1, 1, 1, 1))
F <- c(2, 2)

G <-matrix(ncol = 4, byrow = TRUE,
           data = c(2, 1, 1, 1, -1, 3, 2, 1, -1, 0, 1, 0))
H <- c(-1, 2, 1)
ldei(E, F, G = G, H = H)$X
pars <- lsei(E, F, G = G, H = H, A = diag(nrow = 4), B = rep(0, 4))$X
pars

## -----------------------------------------------------------------------------
(xr <- xranges(E, F, G = G, H = H))

## ----dc, fig.cap = "Fig: Parsimonious solution and ranges of the underdetermined system including equalities and inequalities."----
dotchart(pars, xlim = range(c(pars,xr)), label = paste("x", 1:4, ""))
segments(x0 = xr[ ,1], x1 = xr[ ,2], y0 = 1:nrow(xr), y1 = 1:nrow(xr))

## -----------------------------------------------------------------------------
xs <- xsample(E = E, F = F, G = G, H = H, type = "cda")$X

## ----u2, fig.cap = "Fig: Random sample of the underdetermined system including equalities and inequalities."----
panel.dens <- function(x, ...) {
    USR <- par("usr")
    on.exit(par(usr = USR))
    par(usr = c(USR[1:2], 0, 1.5) )
    DD    <- density(x)
    DD$y  <- DD$y/max(DD$y)
    polygon(DD, col = "grey")
}
xs <- xsample(E = E, F = F, G = G, H = H,  
              jmp = 0.5, sdB = 1)$X
pairs(xs, pch = ".", cex = 2, 
      upper.panel = NULL, diag.panel = panel.dens)

## -----------------------------------------------------------------------------
Va <- c(1, 2, -1, 1)
Vb <- -2
varranges(E, F, G = G, H = H, EqA = Va, EqB = Vb)
summary(varsample(xs, EqA = Va, EqB = Vb))

## -----------------------------------------------------------------------------
A <- matrix(ncol = 4, byrow = TRUE,
            data = c(2, 2, 1, 6, 1, -1, 1, -1))
B <- c(1, 2)

lsei(E, F, G = G, H = H, A = A, B = B)$X

## ----u3, fig.cap = "Fig: Random sample of the underdetermined system including equalities, inequalities, and approximate equations"----
xs <- xsample(E = E, F = F, G = G, H = H, A = A, B = B, 
              jmp = 0.5, sdB = 1)$X
pairs(xs, pch = ".", cex = 2, 
      upper.panel = NULL, diag.panel = panel.dens)

## -----------------------------------------------------------------------------
E <- matrix(ncol = 4, byrow = TRUE,
            data = c(3, 2, 1, 4, 1, 1, 1, 1))
F <- c(2, 2)

G <-matrix(ncol = 4, byrow = TRUE,
           data = c(2, 1, 1, 1, -1, 3, 2, 1, -1, 0, 1, 0))
H <- c(-1, 2, 1)
Cost <- c(1, 2, -1, 4)
linp(E, F, G, H, Cost)

## -----------------------------------------------------------------------------
LP <- linp(E = E, F = F, G = G, H = H, Cost = Cost, ispos = FALSE)
LP$X
LP$solutionNorm

## -----------------------------------------------------------------------------
nn   <- 500

## -----------------------------------------------------------------------------
aa <- runif(nn-1)
bb <- runif(nn)
cc <- runif(nn-1)

## -----------------------------------------------------------------------------
A <-matrix(nrow = nn, ncol = nn, data = 0)
A [cbind((1:(nn-1)),(2:nn))] <- cc
A [cbind((2:nn),(1:(nn-1)))] <- aa
diag(A) <- bb

## -----------------------------------------------------------------------------
abd <- rbind(c(0,cc), bb, c(aa,0))

## -----------------------------------------------------------------------------
A[1:5, 1:5]
aa[1:5]
bb[1:5]
cc[1:5]
abd[ ,1:5]

## -----------------------------------------------------------------------------
B <- runif(nn)
B <- cbind(B, 2*B, 3*B, 4*B, 5*B, 6*B, 7*B, 8*B, 9*B, 10*B)

## -----------------------------------------------------------------------------
print(system.time(
      Full <- solve(A, B) )
     *1000)
print(system.time(
      Band <- Solve.banded(abd, nup = 1, nlow = 1, B))
      *1000)
print(system.time(
      tri  <- Solve.tridiag(aa, bb, cc, B))
      *1000)

## -----------------------------------------------------------------------------
Full[1:5, 1:5]

## -----------------------------------------------------------------------------
head(cbind(Full=Full[,2],Band=Band[,2], Tri=tri[,2]))

## -----------------------------------------------------------------------------
#  0.0  -0.98 -0.79 -0.15                                                  Top
# -1.00  0.25 -0.87  0.35                                                  Top
#  0.78  0.31 -0.85  0.89 -0.69 -0.98 -0.76 -0.82                          blk1
#  0.12 -0.01  0.75  0.32 -1.00 -0.53 -0.83 -0.98
# -0.58  0.04  0.87  0.38 -1.00 -0.21 -0.93 -0.84
# -0.21 -0.91 -0.09 -0.62 -1.99 -1.12 -1.21  0.07
#                          0.78 -0.93 -0.76  0.48 -0.87 -0.14 -1.00 -0.59  blk2
#                         -0.99  0.21 -0.73 -0.48 -0.93 -0.91  0.10 -0.89
#                         -0.68 -0.09 -0.58 -0.21  0.85 -0.39  0.79 -0.71
#                          0.39 -0.99 -0.12 -0.75 -0.68 -0.99  0.50 -0.88
#                                                  0.71 -0.64  0.0   0.48  Bot
#                                                  0.08 100.0 50.00 15.00  Bot

## -----------------------------------------------------------------------------
B <- c(-1.92, -1.27, -2.12, -2.16, -2.27, -6.08, 
       -3.03, -4.62, -1.02, -3.52, 0.55, 165.08)

## -----------------------------------------------------------------------------
Top  <- matrix(nrow = 2, ncol = 4, byrow = TRUE, data =
   c( 0.0,  -0.98, -0.79, -0.15, -1.00,  0.25, -0.87,  0.35))
   
Bot  <- matrix(nrow = 2, ncol = 4, byrow = TRUE, data =
   c( 0.71, -0.64,   0.0,  0.48,  0.08, 100.0, 50.00, 15.00))
   
Blk1 <- matrix(nrow = 4, ncol = 8, byrow = TRUE, data =
   c( 0.78,  0.31, -0.85,  0.89, -0.69, -0.98, -0.76, -0.82,
      0.12, -0.01,  0.75,  0.32, -1.00, -0.53, -0.83, -0.98,
     -0.58,  0.04,  0.87,  0.38, -1.00, -0.21, -0.93, -0.84,
     -0.21, -0.91, -0.09, -0.62, -1.99, -1.12, -1.21,  0.07))
     
Blk2 <- matrix(nrow = 4, ncol = 8, byrow = TRUE, data =
    c( 0.78, -0.93, -0.76,  0.48, -0.87, -0.14, -1.00, -0.59,
      -0.99,  0.21, -0.73, -0.48, -0.93, -0.91,  0.10, -0.89,
      -0.68, -0.09, -0.58, -0.21,  0.85, -0.39,  0.79, -0.71,
       0.39, -0.99, -0.12, -0.75, -0.68, -0.99,  0.50, -0.88))
AR <- array(dim = c(4,8,2), data = c(Blk1, Blk2))

AR

## -----------------------------------------------------------------------------
overlap <- 4
Solve.block(Top, AR, Bot, B, overlap=4)

## -----------------------------------------------------------------------------
B3 <- B
for (i in 2:1000) B3 <- cbind(B3, B*i)

print(system.time(
  X <- Solve.block(Top, AR, Bot, B3, overlap = 4))
  *1000)

X[,c(1, 10, 100, 1000)]

