## ----message = FALSE----------------------------------------------------------
library(nvmix)
library(RColorBrewer)
library(lattice)
doPDF <- FALSE
eval <- TRUE

## ----eval = eval--------------------------------------------------------------
## Generate a random correlation matrix and random limits in dimension d = 5
d <- 5
set.seed(42)
A <- matrix(runif(d * d), ncol = d)
P <- cov2cor(A %*% t(A)) # (randomly generated) correlation matrix
b <-  3 * runif(d) * sqrt(d) # (randomly generated) upper limit
a <- -3 * runif(d) * sqrt(d) # (randomly generated) lower limit

## Specify the mixture distribution parameter
rate <- 1.9 # exponential rate parameter

## Method 1: Use R's qexp() function and provide a list as 'mix'
set.seed(42)
(p1 <- pnvmix(b, lower = a, qmix = list("exp", rate = rate), scale = P))

## Method 2: Define the quantile function manually (note that
##           we do not specify rate in the quantile function here,
##           but conveniently pass it via the ellipsis argument)
set.seed(42)
(p2 <- pnvmix(b, lower = a, qmix = function(u, lambda) -log(1-u)/lambda,
              scale = P, lambda = rate))

## Comparison
stopifnot(all.equal(p1, p2))

## ----eval = eval--------------------------------------------------------------
pnvmix(b, lower = a, qmix = function(u, lambda) -log(1-u)/lambda,
      lambda = rate, scale = P, control = list(pnvmix.abstol = 1e-5))

## ----eval = eval--------------------------------------------------------------
## Define the quantile function of the three-point distribution
## which puts masses 'p' at the numbers 'x'
x <- c(1, 3, 5) # support
p <- c(0.2, 0.3, 0.5) # probabilities
qW <- function(u)
    (u <= p[1]) * x[1] + (u > p[1] & u <= p[1]+p[2]) * x[2] + (u > p[1]+p[2]) * x[3]

## Call pnvmix(); lower defaults to (-Inf,...,-Inf)
set.seed(42)
(p1 <- pnvmix(b, qmix = qW, scale = P))

## ----eval = eval--------------------------------------------------------------
set.seed(42)
p2 <- sum(sapply(1:3, function(k) p[k] * pNorm(b, scale = x[k] * P)))
stopifnot(all.equal(p1, p2, check.attributes = FALSE, tol = 5e-4))

## ----fig.align = "center", fig.width = 7, fig.height = 7, fig.show = "hold", eval = eval----
## Setup
df <- 1.5 # degrees of freedom
maxiter <- 9 # note: i iterations require 3 * 2^8 * 2^i function evaluations
max.fun.evals <- 3 * 2^8 * 2^seq(from = 2, to = maxiter, by = 1)
errors <- matrix(, ncol = length(max.fun.evals), nrow = 4)
nms <- c("Sobol  with preconditioning", "Sobol  w/o  preconditioning",
         "PRNG with preconditioning", "PRNG w/o  preconditioning")
rownames(errors) <- nms

## Computing the errors
## Note:
## - resetting the seed leads to a fairer comparison here
## - set 'verbose' to 0 or FALSE to avoid warnings which inevitably occur
##   due to 'pnvmix.abstol = NULL'
for(i in seq_along(max.fun.evals)) {
    N.max <- max.fun.evals[i]
    ## Sobol with preconditioning
    set.seed(42)
    errors[nms[1],i] <-
        attr(pStudent(b, lower = a, scale = P,  df = df,
                      control = list(pnvmix.abstol = NULL, fun.eval = c(2^6, N.max)),
                      verbose = FALSE), "abs. error")
    ## Sobol without preconditioning
    set.seed(42)
    errors[nms[2],i] <-
        attr(pStudent(b, lower = a, scale = P,  df = df,
                      control = list(precond = FALSE, pnvmix.abstol = NULL,
                                     fun.eval = c(2^6, N.max)),
                      verbose = FALSE), "abs. error")
    ## PRNG with preconditioning
    set.seed(42)
    errors[nms[3],i] <-
        attr(pStudent(b, lower = a, scale = P,  df = df,
                      control = list(method = "PRNG", pnvmix.abstol = NULL,
                                     fun.eval = c(2^6, N.max)),
                      verbose = FALSE), "abs. error")
    ## PRNG without preconditioning
    set.seed(42)
    errors[nms[4],i] <-
        attr(pStudent(b, lower = a, scale = P,  df = df,
                      control = list(method = "PRNG", precond = FALSE,
                                     pnvmix.abstol = NULL, fun.eval = c(2^6, N.max)),
                      verbose = FALSE), "abs. error")
}

## Computing the regression coefficients
coeff <- apply(errors, 1, function(y) lm(log(y) ~ log(max.fun.evals))$coeff[2])
names(coeff) <- nms

## Plot
if(doPDF) pdf(file = (file <- "fig_pnvmix_error_comparison.pdf"),
              width = 7, height = 7)
pal <- colorRampPalette(c("#000000", brewer.pal(8, name = "Dark2")[c(7, 3, 5)]))
cols <- pal(4) # colors
plot(NA, log = "xy", xlim = range(max.fun.evals), ylim = range(errors),
     xlab = "Number of function evaluations", ylab = "Estimated error")
lgnd <- character(4)
for(k in 1:4) {
    lines(max.fun.evals, errors[k,], col = cols[k])
    lgnd[k] <- paste0(nms[k]," (",round(coeff[k], 2),")")
}
legend("topright", bty = "n", lty = rep(1, 4), col = rev(cols), legend = rev(lgnd))
if(doPDF) dev.off()

## ----eval = eval--------------------------------------------------------------
x <- matrix(1:15/15, ncol = d) # evaluation points of the density
set.seed(1)
(d1 <- dnvmix(x, qmix = qW, scale = P)) # computed density values
set.seed(1)
(d2 <- dnvmix(x, qmix = qW, scale = P, log = TRUE)) # log-density
stopifnot(all.equal(d1, exp(d2), check.attributes = FALSE)) # check
## This could have also been obtained via
d3 <- rowSums(sapply(1:3, function(k) p[k] * dNorm(x, scale = x[k] * P)))
stopifnot(all.equal(d1, d3, tol = 1e-10, check.attributes = FALSE))

## ----fig.align = "center", fig.width = 9, fig.height = 6, fig.show = "hold", eval = eval----
n <- 40 # sample size 
x <- matrix(1:n, ncol = 2) # n/2 - two dimensional evaluation points 
m <- mahalanobis(x, center = c(0,0), cov = diag(2)) # for plotting
df <- 2
d3.1 <- dStudent(x, df = df, log = TRUE) # true value
## Specify 'qmix' as function to force estimation of log-density via RQMC
d3.2 <- dnvmix(x, qmix = function(u) 1/qgamma(1-u, shape = df/2, rate = df/2), 
               log = TRUE)
rel.err <- (d3.2 - d3.1) / d3.1
stopifnot(max(abs(rel.err)) < 5e-3) # check 
cols <- pal(2)
if(doPDF) pdf(file = (file <- paste0("fig_dStudentvsdnvmix.pdf")),
              width = 6, height = 6)
plot(sqrt(m), d3.1, type = 'l', col = cols[1], 
     xlab = expression(paste("Mahalanobis Distance ", x^T, x)), ylab = "log-density")
lines(sqrt(m), d3.2, col = cols[2], lty = 2)
legend("topright", c("True log-density", "Estimated log-density"),
       lty = c(1,2), col = cols[1:2], bty = 'n')
if(doPDF) dev.off()

## ----fig.align = "center", fig.width = 6, fig.height = 6, fig.show = "hold", eval = eval----
## Sampling
n <- 500 # sample size
set.seed(42)
r1 <- rnvmix(n, rmix = list("exp", rate = rate)) # uses the default P = diag(2)

## Plot
if(doPDF) pdf(file = (file <- paste0("fig_rnvmix_W_exp.pdf")),
              width = 6, height = 6)
plot(r1, xlab = expression(X[1]), ylab = expression(X[2]))
if(doPDF) dev.off()

## ----fig.align = "center", fig.width = 9, fig.height = 6, fig.show = "hold", eval = eval----
## Sampling
set.seed(42)
r1 <- rnvmix(n, qmix = qW)
r2 <- rnvmix(n, qmix = qW, method = "ghalton")

## Plot
if(doPDF) pdf(file = (file <- paste0("fig_rnvmix_W_three-point.pdf")),
              width = 9, height = 6)
ran <- range(r1, r2)
opar <- par(pty = "s")
layout(t(1:2))
plot(r1, xlab = expression(X[1]), ylab = expression(X[2]),
     main = "Pseudo-random sample", xlim = ran, ylim = ran)
plot(r2, xlab = expression(X[1]), ylab = expression(X[2]),
     main = "Quasi-random sample", xlim = ran, ylim = ran)
layout(1)
par(opar)
if(doPDF) dev.off()

## ----fig.align = "center", fig.width = 6, fig.height = 6, fig.show = "hold", eval = eval----
## Sampling
set.seed(42)
r <- lapply(1:3, function(k) rNorm(p[k] * n, scale = diag(x[k], 2)))

## Plot
if(doPDF) pdf(file = (file <- paste0("fig_rnvmix_W_three-point_via_rNorm.pdf")),
              width = 6, height = 6)
ran  <- range(r)
cols <- pal(4)
opar <- par(pty = "s")
plot(NA, xlim = ran, ylim = ran, xlab = expression(X[1]), ylab = expression(X[2]))
for(k in 1:3) points(r[[k]], col = cols[k+1])
par(opar)
if(doPDF) dev.off() 

## ----fig.align = "center", fig.width = 7, fig.height = 7, fig.show = "hold", eval = eval----
## Sampling
df <- 3.9 # degrees of freedom
factor <- matrix(c(1,0, 0,1, 0,1), ncol = 2, byrow = TRUE) # (3,2)-matrix 'factor'
Sigma <- tcrossprod(factor) # the 'scale' corresponding to factor
stopifnot(Sigma == factor %*% t(factor))
set.seed(42)
r <- rStudent(n, df = df, factor = factor) # sample

## Plot
if(doPDF) pdf(file = (file <- paste0("fig_rnvmix_singular.pdf")),
              width = 6, height = 6)
cloud(r[,3] ~ r[,1] * r[,2], screen = list(z = 115, x = -68),
      xlab = expression(X[1]), ylab = expression(X[2]), zlab = expression(X[3]),
      scales = list(arrows = FALSE, col = "black"),
      par.settings = modifyList(standard.theme(color = FALSE),
                                list(axis.line = list(col = "transparent"),
                                     clip = list(panel = "off"))))
if(doPDF) dev.off()

## ----fig.align = "center", fig.width = 6, fig.height = 6, fig.show = "hold", eval = eval----
set.seed(42) # for reproducibility
## Define 'qmix' as the quantile function of a Par(nu, 1) distribution
qmix <- function(u, nu) (1-u)^(-1/nu)
## Parameters for sampling
n         <- 50
d         <- 3
loc       <- rep(0, d) # true location vector
A         <- matrix(runif(d * d), ncol = d)
scale     <- cov2cor(A %*% t(A)) # true scale matrix
nu        <- 2.4 # true mixing parameter
mix.param.bounds  <- c(1, 10) # nu in [1, 10]
## Sample data using 'rnvmix()':
x <- rnvmix(n, qmix = qmix, nu = nu, loc = loc, scale = scale)
## Call 'fitvnmix()' with 'qmix' as function (so all densities/weights are estimated)
(MyFit21 <- fitnvmix(x, qmix = qmix, mix.param.bounds = mix.param.bounds))
## Call 'fitnvmix()' with 'qmix = "pareto"' in which case an analytical formula
## for the density is used
(MyFit22 <- fitnvmix(x, qmix = "pareto", mix.param.bounds = mix.param.bounds))
stopifnot(all.equal(MyFit21$nu, MyFit22$nu, tol = 5e-2))
## Produce a Q-Q-Plot of the sampled mahalanobis distance versus their theoretical
## quantiles with parameters estimated in 'MyFit21'
if(doPDF) pdf(file = (file <- paste0("fig_fitnvmix_qqplot.pdf")),
              width = 6, height = 6)
qqplot_maha(x, qmix = "pareto", loc = MyFit21$loc, scale = MyFit21$scale,
            alpha = MyFit21$nu)
if(doPDF) dev.off()

