# Copyright (c) Rob Carnell 2026

# More, Garbow, Hillstrom: Testing Unconstrained Optimization Software
# ACM Trans. Math. Software, 7, March 1981, 17--41
# Function as shown in this paper is for optimizing not for non linear equations
# from the paper it is not clear what was used by the authors for their tests
# of non linear equation solvers

# function 27
# Brown almost linear function

brown <- function(x) {
  n <- length(x)
  y <- numeric(n)

  y[1:(n-1)] <- x[1:(n-1)] + sum(x[1:n]) - (n + 1)
  y[n] <- prod(x[1:n]) - 1.0

  y
}

brown_xstart <- function(n) rep(1,n)/2

# function 28
# Discrete boundary value problem

dcbval <- function (x) {
  n <- length(x)
  y <- numeric(n)

  h <- 1/(n+1)

  y <- 2*x + .5*h^2 * (x+1+h*c(1:n))^3
  y[2:n]     <- y[2:n]     - x[1:(n-1)]
  y[1:(n-1)] <- y[1:(n-1)] - x[2:n]

  y
}

dcbval_xstart <- function(n) c(1:n) * (c(1:n)-n-1)/(n+1)^2

################################################################################

# function 29
# Discrete integral equation

dciequ <- function(x) {
  n <- length(x)
  y <- numeric(n)

  h <- 1/(n+1)

  # do it the stupid (and slow?) way

  for (k in 1:n) {
    tk <-  k/(n+1)

    sum1 <- 0
    for (j in 1:k) {
      tj <- j*h
      sum1 <- sum1 + tj * (x[j] + tj + 1.0)^3
    }

    # do this with while since R will count down with for loop
    sum2 <- 0
    j <- k+1
    while(j <= n) {
      tj  <- j*h
      sum2 <- sum2 + (1.0 - tj) * (x[j] + tj + 1.0)^3
      j <- j + 1
    }
    y[k] = x[k] + h*((1.0 - tk) * sum1 + tk * sum2 ) / 2.0
  }

  y
}

dciequ_xstart <- function(n) c(1:n) * (c(1:n)-1-n)/(n+1)^2

################################################################################

# function 7
# Helical valley function

helval <- function(x) {
  c1 <- 0.25
  c2 <- 0.5

  stopifnot(length(x) == 3)
  y <- numeric(3)

  tpi <- 8*atan(1)
  if (x[1] > 0)
    temp1 <- atan(x[2]/x[1])/tpi
  else if (x[1] < 0)
    temp1 <- atan(x[2]/x[1])/tpi + c2
  else
    temp1 <- sign(c1,x[2])

  y[1] <- 10*(x[3] - 10*temp1)
  y[2] <- 10*(sqrt(x[1]^2+x[2]^2) - 1)
  y[3] <- x[3]

  y
}

helval_xstart <- c(-1,0,0)

################################################################################

# Powell singular function
# function 13

pwlsng <- function(x) {
  stopifnot(length(x) == 4)
  y <- numeric(4)
  y[1] <- x[1] + 10*x[2]
  y[2] <- sqrt(5)*(x[3] - x[4])
  y[3] <- (x[2] - 2*x[3])^2
  y[4] <- sqrt(10)*(x[1] - x[4])^2

  y
}

pwlsng_xstart <- c(3, -1, 0, 1)

################################################################################

# Rosenbrock function
# function 1

rosbrk <- function(x) {
  stopifnot(length(x) == 2)
  y <- numeric(2)

  y[1] <- 1 - x[1]
  y[2] <- 10*(x[2] - x[1]^2)

  y
}

rosbrk_xstart <- c(-1.2, 1)

################################################################################

# vardim function
# Function 25
# The function here comes from
#  http://www.zib.de/Numerik/numsoft/CodeLib/codes/newton_testset/problems/VariablyDim.f

# known solutions
vardimsol <- function(n) {
  x <- numeric(n)
  x[1:n] <- 1
  x
}

vardim <- function(x) {
  n <- length(x)
  s <- sum(c(1:n)*(x-1))
  y <- x - 1 + s * (1+2*s^2) * c(1:n)
}

vardiminit <- function (n) {
  x <- (1 - c(1:n) / n)
}

################################################################################

# Watson function
# code adapted from test problems f90 for Minpack
# very difficult function
# function 20

watson <- function(x) {

  n <- length(x)
  y <- numeric(n)

  for(i in 1:29) {
    ti <- i/29

    sum1 <- 0
    t    <- 1
    for (j in 2:n) {
      sum1 <- sum1 + (j-1) * t * x[j]
      t    <- ti * t
    }

    sum2 <- 0
    t    <- 1
    for (j in 1:n) {
      sum2 <- sum2 + t * x[j]
      t    <- ti * t
    }

    t <- 1 / ti
    for(k in 1:n) {
      y[k] <- y[k] + t * (sum1-sum2^2-1) * (k-1-2*ti*sum2)
      t <- t * ti
    }
  }

  y[1] = y[1] + 3*x[1] - 2*x[1]*x[2] + 2*x[1]^3
  y[2] = y[2] + x[2] - x[1]^2 - 1

  y
}

watson_xstart <- function(n) numeric(n)

################################################################################

# Wood function
# function 14

wood <- function(x) {
  c1 <- 2.0e2
  c2 <- 2.02e1
  c3 <- 1.98e1
  c4 <- 1.8e2

  stopifnot(length(x) == 4)
  y <- numeric(4)

  y[1] <- -c1*x[1]*(x[2] - x[1]^2) - (1 - x[1])
  y[2] <- c1*(x[2] - x[1]^2) + c2*(x[2] - 1) + c3*(x[4]-1)
  y[3] <- -c4*x[3]*(x[4] - x[3]^2) - (1 - x[3])
  y[4] <- c4*(x[4] - x[3]^2) + c2*(x[4] - 1) + c3*(x[2]-1)

  y
}

wood_xstart <- c(-3,-1,-3,-1)

################################################################################

# Powell Badly scaled function
# function 3

pwlbsc <- function(x) {
  c1 <- 1.e4
  c2 <- 1.0001

  stopifnot(length(x) == 2)
  y <- numeric(2)
  y[1] <- c1*x[1]*x[2] - 1
  y[2] <- exp(-x[1]) + exp(-x[2]) - c2

  y
}

pwlbscjac <- function(x) {
  stopifnot(length(x) == 2)
  c1 <- 1.e4
  c2 <- 1.0001
  J <- matrix(NA, 2, 2)
  J[1,1] <- c1*x[2]
  J[1,2] <- c1*x[1]
  J[2,1] <- -exp(-x[1])
  J[2,2] <- -exp(-x[2])
  J
}

pwlbsc_xstart <- c(0,1)

################################################################################

# Dennis & Schnabel,1996,"Numerical methods for unconstrained optimization and nonlinear equations", SIAM
# example 6.5.1 page 149

dslnex <- function(x) {
  stopifnot(length(x) == 2)
  y <- numeric(2)
  y[1] <- x[1]^2 + x[2]^2 - 2
  y[2] <- exp(x[1]-1) + x[2]^3 - 2
  y
}

jacdsln <- function(x) {
  stopifnot(length(x) == 2)
  Df <- matrix(NA, 2, 2)
  Df[1,1] <- 2*x[1]
  Df[1,2] <- 2*x[2]
  Df[2,1] <- exp(x[1]-1)
  Df[2,2] <- 3*x[2]^2

  Df
}

# give parameters names which must appear in results
dslnex_xstart <- c(x1=1.5, x2=2)

################################################################################

# Nocedal & Wright (2006) example 11.31

f_nocedal <- function(x) {
  stopifnot(length(x) == 2)
  y <- numeric(2)
  y[1] <- (x[1]+3)*(x[2]^2-7) + 18
  y[2] <- sin(x[2]*exp(x[1])-1)
  y
}

jac_nocedal <- function(x) {
  stopifnot(length(x) == 2)
  J <- matrix(0, nrow=2, ncol=2)
  J[1,1] <- x[2]^2 - 7
  J[1,2] <- (x[1] + 3)*(2*x[2])
  J[2,2] <- cos(x[2]*exp(x[1])-1)*exp(x[1])
  J[2,1] <- x[2]*J[2,2]  # cos(x[2]*exp(x[1])-1)*x[2]*exp(x[1])
  J
}

nocedal_xstart <- c(-1, 2)

################################################################################

# from  BB in Journal of Statistical Software
# from Jstatsoft paper on BB package (Barzila-Borwein)

# @article{Varadhan:Gilbert:2009:JSSOBK:v32i04,
#   author =    "Ravi Varadhan and Paul Gilbert",
#   title = "BB: An R Package for Solving a Large System of Nonlinear Equations and for Optimizing a High-Dimensional Nonlinear Objective Function",
#   journal =   "Journal of Statistical Software",
#   volume =    "32",
#   number =    "4",
#   pages = "1--26",
#   day =   "14",
#   month = "10",
#   year =  "2009",
#   CODEN = "JSSOBK",
#   ISSN =  "1548-7660",
#   bibdate =   "2009-06-29",
#   URL =   "http://www.jstatsoft.org/v32/i04",
#   accepted =  "2009-06-29",
#   acknowledgement = "",
#   keywords =  "",
#   submitted = "2008-12-31",
# }

chandraH <- function(x, c=0.9) {
  n <- length(x)
  k <- 1:n
  mu <- (k - 0.5)/n
  dterm <- outer(mu, mu, function(x1,x2) x1 / (x1 + x2) )
  x - 1 / (1 - c/(2*n) * rowSums(t(t(dterm) * x)))
}

# approximate jacobian
# needed because otherwise it has difficulty in finding a solution

chandraH.jac <- function(x) {
  n <- length(x)
  J <- matrix(0,nrow=n,ncol=n)

  diag(J) <- 1

  J
}

chandra_xstart <- function(n) runif(n)

################################################################################

# last example in
# Hirsch, Pardalos, Resende: Solving systems of nonlinear equations with continuous GRASP
# 0.25 <= x[1] <= 1 and 1.5 <= x[2] <= 2*pi
# two roots

sinexp <- function(x) {
  stopifnot(length(x) == 2)
  f <- numeric(2)
  f[1] <- 0.5*sin(prod(x)) - 0.25*x[2]/pi - 0.5*x[1]
  f[2] <- (1-0.25/pi)*(exp(2*x[1])-exp(1)) + exp(1)*x[2]/pi - 2*exp(1)*x[1]

  f
}

sinexp_xstart <- function(g) {
  xstart <- numeric(2)
  xstart[1] <-  (g*0.25 + 1-g)/2
  xstart[2] <-  (g*1.25 + (1-g)*pi)/2
  xstart
}

################################################################################

# from  BB in Journal of Statistical Software
# from Jstatsoft paper on BB package (Barzila-Borwein)

# @article{Varadhan:Gilbert:2009:JSSOBK:v32i04,
#   author =    "Ravi Varadhan and Paul Gilbert",
#   title = "BB: An R Package for Solving a Large System of Nonlinear Equations and for Optimizing a High-Dimensional Nonlinear Objective Function",
#   journal =   "Journal of Statistical Software",
#   volume =    "32",
#   number =    "4",
#   pages = "1--26",
#   day =   "14",
#   month = "10",
#   year =  "2009",
#   CODEN = "JSSOBK",
#   ISSN =  "1548-7660",
#   bibdate =   "2009-06-29",
#   URL =   "http://www.jstatsoft.org/v32/i04",
#   accepted =  "2009-06-29",
#   acknowledgement = "",
#   keywords =  "",
#   submitted = "2008-12-31",
# }

troesch <- function(x) {
  n <- length(x)
  tnm1 <- 2:(n-1)
  f <- rep(NA, n)

  h <- 1 / (n+1)
  h2 <- 10 * h^2

  f[1]    <- 2 * x[1] + h2 * sinh(10 * x[1]) - x[2]
  f[tnm1] <- 2 * x[tnm1] + h2 * sinh(10 * x[tnm1]) - x[tnm1-1] - x[tnm1+1]
  f[n]    <- 2 * x[n] + h2 * sinh(10 * x[n]) - x[n-1] - 1

  f
}

troesch_xstart <- function(n) runif(n)

################################################################################

# intersection of circles
# http://paulbourke.net/geometry/circlesphere/

#    The following note describes how to find the intersection point(s) between two circles on a plane,
#    the following notation is used. The aim is to find the two points P3 = (x3, y3) if they exist.
#
#    see circle_2circle1.gif

#    First calculate the distance d between the center of the circles. d = ||P1 - P0||.
#
#    If d > r0 + r1 then there are no solutions, the circles are separate.
#    If d < |r0 - r1| then there are no solutions because one circle is contained within the other.
#    If d = 0 and r0 = r1 then the circles are coincident and there are an infinite number of solutions.
#    Considering the two triangles P0P2P3 and P1P2P3 we can write
#    a2 + h2 = r02 and b2 + h2 = r12
#    Using d = a + b we can solve for a,
#
#    a = (r02 - r12 + d2 ) / (2 d)
#    It can be readily shown that this reduces to r0 when the two circles touch at one point, ie: d = r0 + r1
#
#    Solve for h by substituting a into the first equation, h2 = r02 - a2
#
#    So
#    P2 = P0 + a ( P1 - P0 ) / d
#    And finally, P3 = (x3,y3) in terms of P0 = (x0,y0), P1 = (x1,y1) and P2 = (x2,y2), is
#
#    x3 = x2 +- h ( y1 - y0 ) / d
#    y3 = y2 -+ h ( x1 - x0 ) / d


# http://www.ambrsoft.com/TrigoCalc/Circles2/Circle2.htm

# (x-a)^2 + (y-b)^2 = r0^2
# (x-c)^2 + (y-d)^2 = r1^2

# D = sqrt((c-a)^2+(d-b)^2)

# r0+r1> D and D> abs(r0-r1)

# line connecting two intersection points

# y =(a-c)/(d-b) * x +( (r0^2-r1^2) + (c^2-a^2) + (d^2-b^2) ) / ( 2*(d-b) )

# equation of line connecting two ceneters

# y = (d-b)/(c-a) * x - a*(d-b)/(c-a) + b

# Delta=0.25 * sqrt( (D+r0+r1)*(D+r0-r1)*(D-r0+r1)*(-D+r0+r1) )

# x1,2 = (a+c)/2 + (c-a)*(r0^2-r1^2)/(2*D^2) +/- 2*(b-d)/D^2 * Delta
# y1,2 = (b+d)/2 + (d-b)*(r0^2-r1^2)/(2*D^2) -/+ 2*(a-c)/D^2 * Delta

# (x-1)^2 + (y-2)^2 = 9
# (x-3)^2 + (y+1)^2 = 16
# a=1, b=2,r0=3
# c=3,d=-1,r1=4

# intersection points (3.86,2.910) and (-0.94,-0.29)

xc1 <- c(1,2)
xc2 <- c(3,-1)
r0 <- 3
r1 <- 4

circle.intersect <- function(xc1,xc2,r0,r1) {
  a <- xc1[1]
  b <- xc1[2]
  cz <- xc2[1]
  d <- xc2[2]

  D <- sqrt((cz-a)^2+(d-b)^2)
  Delta <- 0.25 * sqrt( (D+r0+r1)*(D+r0-r1)*(D-r0+r1)*(-D+r0+r1) )

  x1 <- (a+cz)/2 + (cz-a)*(r0^2-r1^2)/(2*D^2) + 2*(b-d)/D^2 * Delta
  x2 <- (a+cz)/2 + (cz-a)*(r0^2-r1^2)/(2*D^2) - 2*(b-d)/D^2 * Delta

  y1 <- (b+d)/2 + (d-b)*(r0^2-r1^2)/(2*D^2) - 2*(a-cz)/D^2 * Delta
  y2 <- (b+d)/2 + (d-b)*(r0^2-r1^2)/(2*D^2) + 2*(a-cz)/D^2 * Delta

  return(list(P1=c(x1,y1),P2=c(x2,y2)))
}

fcircle <- function(x,xc1,xc2,r0,r1) {
  f <- numeric(2)
  f[1] <- (x[1]-xc1[1])^2+(x[2]-xc1[2])^2 - r0^2
  f[2] <- (x[1]-xc2[1])^2+(x[2]-xc2[2])^2 - r1^2
  return(f)
}

################################################################################

# Steady-State solution for reaction rate equations
# Shacham homotopy method (discrete changing of one or more parameters)
# M. Shacham: Numerical Solution of Constrained Non-linear algebriac equations
# International Journal for Numerical Methods in Engineering, 1986, pp.1455--1481.

# solution should always be > 0
# second initial point converges to infeasible solution

# Problem 2, page 1463/1464

cutlip <- function(x, k1, k2, k3, kr1, kr2) {
  stopifnot(length(x) == 6)

  r <- numeric(6)

  r[1] = 1 - x[1] - k1*x[1]*x[6] + kr1 * x[4]
  r[2] = 1 - x[2] - k2*x[2]*x[6] + kr2 * x[5]
  r[3] = -x[3] + 2*k3*x[4]*x[5]
  r[4] = k1*x[1]*x[6] - kr1*x[4] - k3*x[4]*x[5]
  r[5] = 1.5*(k2*x[2]*x[6] - kr2*x[5]) - k3*x[4]*x[5]
  r[6] = 1 - x[4] - x[5] - x[6]

  r
}

cutlip_xstart <- function(n) runif(n, 0, 1)

################################################################################

# R. Baker Kearfott, Some tests of Generalized Bisection,
# ACM Transactions on Methematical Software, Vol. 13, No. 3, 1987, pp 197-220

# Robot kinematics example Kearfott (4.2 Problem 11)

a <- c(4.371e-3, -0.3578   , -0.1238,
       -1.637e-3, -0.9338   ,  1.0,
       -0.3571  ,  0.2238   ,  0.7623,
       0.2638  , -0.7745e-1, -0.6734,
       -0.6022  ,  1.0      ,  0.3578,
       4.731e-3, -0.7623   ,  0.2238,
       0.3461
)

kearfott <- function(x) {
  stopifnot(length(x) == 8)
  y <- numeric(8)
  y[1] <- a[1]*x[1]*x[3]  + a[2]*x[2]*x[3] + a[3]*x[1] + a[4]*x[2]  + a[5]*x[4]  + a[6]*x[7] + a[7]
  y[2] <- a[8]*x[1]*x[3]  + a[9]*x[2]*x[3] + a[1]*x[1] + a[11]*x[2] + a[12]*x[4] + a[13]
  y[3] <- a[14]*x[6]*x[8] + a[15]*x[1]     + a[16]*x[2]
  y[4] <- a[17]*x[1] + a[18]*x[2] + a[19]
  y[5] <- x[1]^2 + x[2]^2-1
  y[6] <- x[3]^2 + x[4]^2-1
  y[7] <- x[5]^2 + x[6]^2-1
  y[8] <- x[7]^2 + x[8]^2-1
  y
}

kearfott_xstart <- runif(8, -1, 1)

################################################################################

# A high-degree polynomial system (section 4.3 Problem 12)
# There are 12 real roots (and 126 complex roots to this system!)

hdp <- function(x) {
  f <- numeric(length(x))
  f[1] <- 5 * x[1]^9 - 6 * x[1]^5 * x[2]^2 + x[1] * x[2]^4 + 2 * x[1] * x[3]
  f[2] <- -2 * x[1]^6 * x[2] + 2 * x[1]^2 * x[2]^3 + 2 * x[2] * x[3]
  f[3] <- x[1]^2 + x[2]^2 - 0.265625
  f
}

hdp_xstart <- runif(3, -1, 1)

# Two intersecting parabolas (4.3 Problem 14)

twoip <- function(x) {
  stopifnot(length(x) == 2)
  f <- numeric(length(x))
  f[1] <- x[1]^2 - 4*x[2]
  f[2] <- x[2]^2 - 2*x[1] + 4*x[2]
  f
}

twoip_xstart <- runif(2, -1, 1)

################################################################################

# http://math.fullerton.edu/mathews/n2003/BroydenMethodMod.html
# page can only be found by searching for broyden site:*.fullerton.edu

# Exercize 1 (2 solutions)

broyden_f1 <- function(par) {
  stopifnot(length(par) == 2)
  x <- par[1]
  y <- par[2]
  f <- numeric(2)

  f[1] <- 1 - 4*x + 2*x^2 - 2*y^3
  f[2] <- -4 + x^4 + 4*y + 4*y^4

  f
}

broyden_f2 <- function(par) {
  stopifnot(length(par) == 3)
  x <- par[1]
  y <- par[2]
  z <- par[3]
  f <- numeric(3)

  f[1] <- 9*x^2 + 36*y^2 + 4*z^2 - 36
  f[2] <- x^2 - 2*y^2 - 20*z
  f[3] <- x^2 - y^2 + z^2

  f
}

################################################################################

# Exercise 7 (1 solution)
broyden_f7 <- function(par) {
  stopifnot(length(par) == 4)
  w <- par[1]
  x <- par[2]
  y <- par[3]
  z <- par[4]
  f <- numeric(4)

  f[1] <- w^2 + x^2 + 3*y^2 - z^3 - 5
  f[2] <- w + x^3 - 2*y^2 - 10*z
  f[3] <- 20 - w + x^2 + y^3 + z^2
  f[4] <- w^3 + x - y^3 + z - 10

  f
}

broyden_f7jac <- function(par) {
  stopifnot(length(par) == 4)
  w <- par[1]
  x <- par[2]
  y <- par[3]
  z <- par[4]
  J <- matrix(NA,nrow=4,ncol=4)
  J[1,1] <- 2*w   ;  J[1,2] <- 2*x  ;    J[1,3] <- 6*y   ;    J[1,4] <- -3*z^2
  J[2,1] <- 1     ;  J[2,2] <- 3*x^2;    J[2,3] <- -4*y  ;    J[2,4] <- -10
  J[3,1] <- -1    ;  J[3,2] <- 2*x  ;    J[3,3] <- 3*y^2 ;    J[3,4] <- 2*z
  J[4,1] <- 3*w^2 ;  J[4,2] <- 1    ;    J[4,3] <- -3*y^2;    J[4,4] <- 1

  J
}

broyden_f7_xstart <- c(-4,4,-4,4)

################################################################################

# Minoru Kuno and J.D. Seader:
# Computing All Real Solutions to Systems of Nonlinear Equations with a Global Fixed-Point Homotopy
# Industrial & Engineering Chemistry Research,  Vol 1988, 27, pp. 1320-1329

# test Kuno Seader equation set 17

kunoseader_f <- function(x) {
  stopifnot(length(x) == 2)
  y <- numeric(2)
  y[1] <- 2*x[1]^3 + 2 * x[1]*x[2] - 22*x[1] + x[2]^2 + 13
  y[2] <- x[1]^2 + 2*x[1]*x[2]+2*x[2]^3 - 14 * x[2] + 9
  y
}

# test Kuno Seader equation set 34

kunoseader_f34 <- function(x) {
  stopifnot(length(x) == 3)
  y <- numeric(3)
  y[1] <- (x[1]-x[2]^2)*(x[1]-sin(x[2]))
  y[2] <- (cos(x[2])-x[1])*(x[2]-cos(x[1]))
  # x[3] slack variable to ensure 0 <= x[1]/[2] <= 1
  y[3] <- x[2]*(x[2]-1)+x[3]^2
  y
}

################################################################################

# Becker Rustem function

# Robin Becker, Berc Rustem: Algorithms for solving nonlinear dynamic decision models,
#                            Annals of Operations Research 44 , 1993, pp. 117-42
# example 1A (function, start)

bkrust <- function(x) {
  n <- length(x)
  y <- numeric(n)

  y[1] <- 4*(x[1]-x[2]^2)

  k <- 2:(n-1)
  y[k] <- 8*x[k]*(x[k]^2-x[k-1]) - 2*(1-x[k]) + 4*(x[k]-x[k+1]^2)

  y[n] <- 8*x[n]*(x[n]^2 - x[n-1]) - 2 * (1 - x[n])
  y
}

bkrjac <- function(x) {
  n <- length(x)
  J <- matrix(numeric(n*n),n,n)

  J[1,1] <- 4
  J[1,2] <- -8*x[2]
  J[n,n-1] <- -8*x[n]
  J[n,n  ] <- 24*x[n]^2 - 8*x[n-1]+2
  for (p in 2:(n-1)) {
    J[p,p-1] <- -8*x[p]
    J[p,p]   <- 24*x[p]^2 - 8*x[p-1] + 6
    J[p,p+1] <- -8*x[p+1]
  }
  J
}

bkrust_xstart <- rep(-2, 36)

################################################################################

# Boggs function
# mentioned in http://folk.uib.no/ssu029/Pdf_file/Testproblems/testprobRheinboldt03.pdf
# problem 6
# mentioned in Moon Sup Song, GENERALIZED ADAMS METHODS FOR THE SOLUTION OF SYSTEMS OF NONLINEAR EQUATIONS,
#              Bulletin Korean Mathematical Society, Vol, 14, No. 1, 1977

# Rheinboldt mentions: 3 solutions: c(-1,2), c(1,0), c(0.5sqrt(2),2)

boggs <- function(x, a=0.5*pi) {
  y <- c(0,0)
  y[1] <- x[1]^2 - x[2] + 1.0
  y[2] <- x[1] - cos(a * x[2])
  return(y)
}

boggs_xstart <- c(1,0)

################################################################################

# Economic modelling problem
#  A. P. Morgan, Solving Polynomial Systems Using Continuation for Scientiﬁc
#                and Engineering Problems. Englewood Cliffs, NJ: Prentice-Hall, 1987.

econ <- function(x, const) {
  n <- length(x)
  y <- numeric(n)

  for( k in 1:(n-1) ) {
    y[k] <- (x[k] + sum(x[1:(n-k-1)]*x[(k+1):(n-1)]))*x[n] - const[k]
  }
  y[n] <- sum(x[1:(n-1)]) + 1
  y
}

################################################################################

# Hiebert 3rd chemical reaction function;
# used in
# M. Shacham: Numerical Solution of Constrained Non-linear algebriac equations
# International Journal for Numerical Methods in Engineering, 1986, pp.1455--1481.

hiebert <- function(x, R=10) {
  sqrtz <- function(x) { tryCatch(sqrt(x), warning=function(w) Inf, error=function(e) Inf) }

  stopifnot(length(x) == 10)

  y <- numeric(length(x))
  xs <- sum(x)

  y[1] <- x[1]+x[4]-3
  y[2] <- 2*x[1]+x[2]+x[4]+x[7]+x[8]+x[9]+2*x[10] - R
  y[3] <- 2*x[2]+2*x[5]+x[6]+x[7]-8
  y[4] <- 2*x[3]+x[5]-4*R
  y[5] <- x[1]*x[5]-0.193*x[2]*x[4]
  # y[6] <- x[6]*x[2]^0.5 - 0.002597*(x[2]*x[4]*xs)^0.5
  # y[7] <- x[7]*x[4]^0.5 - 0.003448*(x[1]*x[4]*xs)^0.5
  y[6] <- x[6]*sqrtz(x[2]) - 0.002597*sqrtz(x[2]*x[4]*xs)
  y[7] <- x[7]*sqrtz(x[4]) - 0.003448*sqrtz(x[1]*x[4]*xs)
  y[8] <- x[8]*x[4] - 1.799e-5*x[2]*xs
  y[9] <- x[9]*x[4] - 2.155e-4*x[1]*sqrtz(x[3]*xs)
  y[10]<- x[10]*x[4]^2 - 3.846e-5*x[4]^2*xs

  y
}

################################################################################

# li five-diagonal function

# Example 6.3 from  Guangye Li, Succesive column correction algorithms for solving sparse nonlinear systems of equations,
#                   Mathematical Programming 43, 1989,, pp. 187-207

li5 <- function(x) {
  n <- length(x)
  y <- numeric(n)

  y[1] <- 4*(x[1]-x[2]^2) + x[2] - x[3]^2
  y[2] <- 8*x[2]*(x[2]^2-x[1]) -2*(1-x[2])+4*(x[2]-x[3]^2)+x[3]-x[4]^2

  k <- 3:(n-2)
  y[k] <- 8*x[k]*(x[k]^2-x[k-1]) - 2*(1-x[k]) + 4*(x[k]-x[k+1]^2)+x[k-1]^2-x[k-2]+x[k+1]-x[k+2]^2

  y[n-1] <-  8*x[n-1]*(x[n-1]^2-x[n-2])-2*(1-x[n-1])+4*(x[n-1]-x[n]^2) +x[n-2]^2-x[n-3]
  y[n] <- 8*x[n]*(x[n]^2 - x[n-1]) - 2 * (1 - x[n]) + x[n-1]^2-x[n-2]
  y
}

li5jac <-  function(x) {
  n <- length(x)
  J <- matrix(numeric(n*n),n,n)

  J[1,1] <- 4
  J[1,2] <- -8*x[2] + 1
  J[1,3] <- -2*x[3]

  J[2,1] <- -8*x[2]
  J[2,2] <- 24*x[2]^2 -8*x[1] + 6
  J[2,3] <- -8*x[3] + 1
  J[2,4] <- -2*x[4]

  for (p in 3:(n-2)) {
    J[p,p-2] <- -1
    J[p,p-1] <- -8*x[p] + 2*x[p-1]
    J[p,p]   <- 24*x[p]^2 - 8*x[p-1] + 6
    J[p,p+1] <- -8*x[p+1]  +1
    J[p,p+2] <- -2*x[p+2]
  }

  J[n-1,n-3]  <- -1
  J[n-1,n-2]  <- -8*x[n-1] + 2*x[n-2]
  J[n-1,n-1]  <-  24*x[n-1]^2 - 8*x[n-2] + 6
  J[n-1,n  ]  <- -8*x[n]
  J[n,n-2] <-  -1
  J[n,n-1] <- -8*x[n] + 2*x[n-1]
  J[n,n  ] <- 24*x[n]^2 - 8*x[n-1]+2

  J
}

li5_xstart <- rep(-2, 36)

################################################################################

# examples of systems of nonlinear equations from
# A. Neumaier: Introduction to Numerical Analysis, Cambridge University Press, 2001
# Example 6.2.3. (page 314)

neum1 <- function(x) {

  y = numeric(length(x))
  y[1] = x[2]^2 - 3*x[1]^2
  y[2] = x[1]^2+x[1]*x[3]+x[3]^2 - 3*x[2]^2
  y[3] = x[2]^2 + x[2] + 1 - 3*x[3]^2

  y
}

neum1jac <- function(x) {

  n <- length(x)
  J <- matrix(numeric(n*n),n,n)

  J[1,1] = -6*x[1]
  J[1,2] =  2*x[2]
  J[2,1] = 2*x[1] + x[3]
  J[2,2] = -6*x[2]
  J[2,3] = x[1]+2*x[3]
  J[3,2] = 2*x[2]+1
  J[3,3] = -6*x[3]

  J
}

neum1_xstart <- c(0.25,0.50,0.75)

################################################################################

# Example system from Octave manual

octxmpl <- function(x) {
  y <- numeric(2)
  y[1] <- -2*x[1]^2 + 3 * x[1]*x[2]   + 4 * sin(x[2]) - 6
  y[2] <-  3*x[1]^2 - 2 * x[1]*x[2]^2 + 3 * cos(x[1]) + 4

  y
}

octxmpl_xstart <- c(1, 2)

################################################################################

#    Example 2 from
#     W.C. Rheinboldt and J.V.Burkardt
#        A Locally Parameterized Continuation Process
#        AMS Transactions on Mathematical Software
#        June 1983 Vol.9 No. 2

#    Maneuvering airplanes at high altitude
#     x1 roll  rate
#     x2 pitch rate
#     x3 yaw   rate
#     x4 angle of attack
#     x5 side slip angle
#     x6 elevator angle
#     x7 aileron angle
#     x8 rudder angle

# very very difficult
# you need a maximum stepsize to really solve this properly
# not available here

pcex <- function(x) {
  stopifnot(length(x) == 5)
  y <- numeric(5)

  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  x4 <- x[4]
  x5 <- x[5]
  x6 <- 0.0
  x7 <- -0.18691
  x8 <- 0.0

  phi1 <- -0.727 * x2 * x3 + 8.39 * x3 * x4 - 684.4 * x4*x5 + 63.5 * x4 * x7

  phi2 <-  0.949 * x1 * x3 + 0.173 * x1 * x5
  phi3 <- -0.716 * x1 * x2 - 1.578 * x1 * x4 + 1.132 * x4*x7
  phi4 <- - x1 * x5
  phi5 <-   x1 * x4

  a1 <- -3.933 * x1 + 0.107 * x2 + 0.126 * x3 - 9.99 * x5 - 45.83 * x7 -7.64 * x8
  a2 <- -0.987 * x2 - 22.95 * x4 - 28.37 * x6
  a3 <-  0.002 * x1 - 0.235 * x3 + 5.67 * x5 - 0.921 * x7 - 6.51 * x8
  a4 <-  x2 - x4 - 0.168 * x6
  a5 <- -x3 - 0.196 * x5 - 0.0071 * x7

  y[1] <- a1 + phi1
  y[2] <- a2 + phi2
  y[3] <- a3 + phi3
  y[4] <- a4 + phi4
  y[5] <- a5 + phi5

  y
}

pcexjac <- function(x) {
  stopifnot(length(x) == 5)
  n  <- length(x)
  df <- matrix(numeric(n*n),n,n)

  x1 <- x[1]
  x2 <- x[2]
  x3 <- x[3]
  x4 <- x[4]
  x5 <- x[5]
  x6 <- 0.0
  x7 <- -0.18691
  x8 <- 0.0

  df[1,1] <- -3.933
  df[1,2] <-  -0.727 * x3 + 0.107
  df[1,3] <-  -0.727 * x2 + 8.39 * x4 + 0.126
  df[1,4] <-  8.39 * x3 - 684.4 * x5 + 63.5 * x7
  df[1,5] <- -684.4 * x4 - 9.99

  df[2,1] <- 0.949 * x3 + 0.173 * x5
  df[2,2] <- -0.987
  df[2,3] <- 0.949 * x1
  df[2,4] <- -22.95
  df[2,5] <- 0.173 * x1

  df[3,1] <- -0.716 * x2 - 1.578 * x4 + 0.002
  df[3,2] <- -0.716 * x1
  df[3,3] <- -0.235
  df[3,4] <- -1.578 * x1 + 1.132 * x7
  df[3,5] <- 5.67

  df[4,1] <- -x5
  df[4,2] <- 1.0
  df[4,3] <- 0.0
  df[4,4] <- -1.0
  df[4,5] <- -x1

  df[5,1] <- x4
  df[5,2] <- 0.0
  df[5,3] <- -1.0
  df[5,4] <- x1
  df[5,5] <- -0.196

  df
}

pcex_xstart <- c(2.5, -0.2, 0.05, 0.01, -0.1 )

pcex_xsol <- c( 2.587674,  -0.223680,   0.054719,   0.013685,  -0.091730 )

################################################################################

# Fractional conversion in a chemical reactor
# Shacham homotopy method (discrete changing of one or more parameters)
# M. Shacham: Numerical Solution of Constrained Non-linear algebriac equations
# International Journal for Numerical Methods in Engineering, 1986, pp.1455--1481.

# solution should always be > 0

# Problem 3, page 1466

# 0 <= x(1) <=1
# function is undefined for 0.8<=x<=1
# radius of convergence is very small: 0.705<=x<=0.799

frac1 <- function(x)   {
    y <- numeric(1)
	y[1] <- x[1]/(1-x[1])-5*log(0.4*(1-x[1])/(0.4-0.5*x[1])) + 4.45977
	y
}

# transformed equation
# all x[i] >= 0

frac2 <- function(x){
    y <- numeric(3)
	y[1] = x[1] / x[2] - 5 * log(0.4*x[2]/x[3]) + 4.45977
	y[2] = x[2] + x[1] - 1
	y[3] = x[3] + 0.5*x[1] - 0.4
	y
}

