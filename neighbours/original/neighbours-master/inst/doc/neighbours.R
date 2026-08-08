### R code from vignette source 'neighbours.Rnw'

###################################################
### code chunk number 1: options
###################################################
## store current options and graphics settings,
## to be restored at the end of the vignette
old.options <- options(continue = "  ", digits = 3,
                       width = 60, useFancyQuotes = FALSE)
old.par <- par(no.readonly = TRUE)
pv <- packageVersion("neighbours")
pv <- gsub("(.*)[.](.*)", "\\1-\\2", pv)


###################################################
### code chunk number 2: package-seed
###################################################
library("neighbours")
set.seed(3477)


###################################################
### code chunk number 3: LSopt
###################################################
LSopt <- if (requireNamespace("NMOF")) {
             NMOF::LSopt
         } else
             function(OF, algo = list(), ...) {
                 xc  <- algo$x0
                 xcF <- OF(xc, ...)
                 for (s in seq_len(algo$nI)) {
                     xn <- algo$neighbour(xc, ...)
                     xnF <- OF(xn, ...)
                     if (xnF <= xcF) {
                         xc  <- xn
                         xcF <- xnF
                     }
                 }
                 list(xbest = xc, OFvalue = xcF)
             }


###################################################
### code chunk number 4: data
###################################################
ny <- 50     ## length of y, number of rows of X
nx <- 500    ## number of columns of X
y <- runif(ny)
X <- array(runif(ny * nx), dim = c(ny, nx))


###################################################
### code chunk number 5: x0
###################################################
x0 <- logical(nx)
x0[1:15] <- TRUE
head(x0, 20)


###################################################
### code chunk number 6: objective-fun
###################################################
column_cor <- function(x, X, y)
    cor(rowMeans(X[, x]), y)


###################################################
### code chunk number 7: cor-x0
###################################################
column_cor(x0, X, y)


###################################################
### code chunk number 8: nb
###################################################
nb <- neighbourfun(type = "logical", kmin = 10, kmax = 20)


###################################################
### code chunk number 9: LSopt-run
###################################################
sol.ls <- LSopt(column_cor,
                list(neighbour = nb,
                     x0 = x0,      ## initial solution
                     nI = 3000,    ## number of iterations
                     printBar = FALSE),
                X = X, y = y)


###################################################
### code chunk number 10: check-sol
###################################################
column_cor(sol.ls$xbest, X, y)


###################################################
### code chunk number 11: constraints-check
###################################################
sum(sol.ls$xbest)


###################################################
### code chunk number 12: neighbours.Rnw:170-186
###################################################
par(mfrow = c(1, 2), las = 1, bty = "n",
    mar = c(3, 3, 1, 0.5), mgp = c(1.75, 0.25, 0),
    tck = 0.02, cex = 0.7)
plot(y, rowMeans(X[, x0]),
     main = "Initial solution",
     pch = 19, cex = 0.5,
     ylim = c(0.3, 0.7),
     xlab = "y",
     ylab = "Mean of linear combination of columns")
par(yaxt = "n")
plot(y, rowMeans(X[, sol.ls$xbest]),
     main = "Result of Local Search",
     pch = 19, cex = 0.5,
     ylim = c(0.3, 0.7),
     xlab = "y")
axis(4)


###################################################
### code chunk number 13: neighbours.Rnw:211-214
###################################################
x <- logical(9L)
x[4:6] <- TRUE
compare_vectors(x)


###################################################
### code chunk number 14: restrict
###################################################
active <- c(rep(FALSE, 3),
            rep(TRUE, length(x) - 3))
active
nb <- neighbourfun(type = "logical", kmin = 3, kmax = 3, active = active)


###################################################
### code chunk number 15: neighbours.Rnw:231-237
###################################################
xs <- list()
xs[[1]] <- x
for (i in 1:10) {
    xs[[length(xs) + 1]] <- x <- nb(x)
}
do.call(compare_vectors, xs)


###################################################
### code chunk number 16: variance1
###################################################
variance <- function(x, S, ...)
    x %*% S %*% x


###################################################
### code chunk number 17: variance2
###################################################
variance2 <- function(x, R, ...)
    var(R %*% x)


###################################################
### code chunk number 18: create-R
###################################################
if (!requireNamespace("NMOF")) {
    R <- array(rnorm(120*20, sd = 0.03), dim = c(120, 20))
} else
    R <- NMOF::randomReturns(na = 20, 120, 0.03, rho = 0.6)
S <- cov(R)  ## Sigma
x0 <- rep(1/ncol(R), ncol(R))


###################################################
### code chunk number 19: nb
###################################################
nb <- neighbourfun(type = "numeric",
                   min = 0, max = 0.2,
                   stepsize = 0.005)


###################################################
### code chunk number 20: mv-ls1
###################################################
sol.ls <- LSopt(variance,
                list(x0 = x0,
                     neighbour = nb,
                     nI = 2000,
                     printBar = FALSE),
                S = S)

sol.qp <- if (requireNamespace("NMOF") &&
              requireNamespace("quadprog"))
              round(NMOF::minvar(S, wmin = 0, wmax = 0.2), 2) else NA
data.frame(LS = round(sol.ls$xbest, 2)[1:10],
           QP = sol.qp[1:10])


###################################################
### code chunk number 21: mv-ls2
###################################################
sol.ls2 <- LSopt(variance2,
                list(x0 = x0,
                     neighbour = nb,
                     nI = 2000,
                     printBar = FALSE),
                R = R)

data.frame(LS2 = round(sol.ls2$xbest, 2)[1:10],
           QP = sol.qp[1:10])


###################################################
### code chunk number 22: neighbours.Rnw:358-363
###################################################
semivariance <- function(x, R, ...) {
    Rx <- R %*% x
    Rx.lower <- pmin(Rx - mean(Rx), 0)
    sum(Rx.lower)
}


###################################################
### code chunk number 23: attr
###################################################
attr(x0, "Ax") <- R %*% x0


###################################################
### code chunk number 24: variance3
###################################################
variance3 <- function(x, ...)
    var(attr(x, "Ax"))


###################################################
### code chunk number 25: nb-update
###################################################
nb_upd <- neighbourfun(type = "numeric",
                       min = 0, max = 0.2,
                       stepsize = 0.005,
                       update = "Ax", A = R)


###################################################
### code chunk number 26: mv-ls3
###################################################
sol.ls3 <- LSopt(variance3,
                list(x0 = x0,
                     neighbour = nb_upd,
                     nI = 2000,
                     printBar = FALSE))

data.frame(LS = round(sol.ls$xbest,  2)[1:10],
           LS2 = round(sol.ls2$xbest, 2)[1:10],
           LS3 = round(sol.ls3$xbest, 2)[1:10],
           QP = sol.qp[1:10])


###################################################
### code chunk number 27: restore-options
###################################################
## restore options and graphics settings
par(old.par)
options(old.options)


