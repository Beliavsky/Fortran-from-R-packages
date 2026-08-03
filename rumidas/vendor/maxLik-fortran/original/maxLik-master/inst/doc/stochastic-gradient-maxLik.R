### R code from vignette source 'stochastic-gradient-maxLik.Rnw'

###################################################
### code chunk number 1: foo
###################################################
options(keep.source = TRUE, width = 60,
        try.outFile=stdout()  # make try to produce error messages
        )
foo <- packageDescription("maxLik")


###################################################
### code chunk number 2: stochastic-gradient-maxLik.Rnw:163-164 (eval = FALSE)
###################################################
## maxSGA(fn, grad, start, nObs, control)


###################################################
### code chunk number 3: stochastic-gradient-maxLik.Rnw:375-379
###################################################
i <- which(names(MASS::Boston) == "medv")
X <- as.matrix(MASS::Boston[,-i])
X <- cbind("const"=1, X)  # add constant
y <- MASS::Boston[,i]


###################################################
### code chunk number 4: stochastic-gradient-maxLik.Rnw:384-385
###################################################
colMeans(X)


###################################################
### code chunk number 5: stochastic-gradient-maxLik.Rnw:392-393
###################################################
eigenvals <- eigen(crossprod(X))$values


###################################################
### code chunk number 6: stochastic-gradient-maxLik.Rnw:404-407
###################################################
betaX <- solve(crossprod(X)) %*% crossprod(X, y)
betaX <- drop(betaX)  # matrix to vector
betaX


###################################################
### code chunk number 7: stochastic-gradient-maxLik.Rnw:431-436
###################################################
gradloss <- function(theta, index)  {
   e <- y[index] - X[index,,drop=FALSE] %*% theta
   g <- t(e) %*% X[index,,drop=FALSE]
   2*g/length(index)
}


###################################################
### code chunk number 8: gradonly
###################################################
library(maxLik)
set.seed(3)
start <- setNames(rnorm(ncol(X), sd=0.1), colnames(X))
                           # add names for better reference
res <- try(maxSGA(grad=gradloss,
           start=start,
           nObs=nrow(X),
           control=list(iterlim=1000)
           )
    )


###################################################
### code chunk number 9: stochastic-gradient-maxLik.Rnw:476-483
###################################################
res <- maxSGA(grad=gradloss,
              start=start,
              nObs=nrow(X),
              control=list(iterlim=1000,
                           SG_clip=1e4)  # limit ||g|| <= 100
              )
summary(res)


###################################################
### code chunk number 10: stochastic-gradient-maxLik.Rnw:495-499
###################################################
loss <- function(theta, index) {
   e <- y[index] - X[index,] %*% theta
   -crossprod(e)/length(index)
}


###################################################
### code chunk number 11: stochastic-gradient-maxLik.Rnw:514-538
###################################################
res <- maxSGA(loss, gradloss,
              start=start,
              nObs=nrow(X),
              control=list(iterlim=1000,
                           # will misbehave with larger numbers
                           SG_clip=1e4,
                           SG_learningRate=0.001,
                           storeParameters=TRUE,
                           storeValues=TRUE
                           )  
              )
par <- storedParameters(res)
val <- storedValues(res)
par(mfrow=c(1,2))
plot(par[,1], par[,2], type="b", pch=".",
     xlab=names(start)[1], ylab=names(start)[2], main="Parameters")
## add some arrows to see which way the parameters move
iB <- c(40, nrow(par)/2, nrow(par))
iA <- iB - 10
arrows(par[iA,1], par[iA,2], par[iB,1], par[iB,2], length=0.1)
##
plot(seq(length=length(val))-1, -val, type="l",
     xlab="epoch", ylab="MSE", main="Loss",
     log="y")


###################################################
### code chunk number 12: stochastic-gradient-maxLik.Rnw:565-570
###################################################
i <- sample(nrow(X), 0.8*nrow(X))  # training indices, 80% of data
Xt <- X[i,]  # training data
yt <- y[i]
Xv <- X[-i,]  # validation data
yv <- y[-i]


###################################################
### code chunk number 13: stochastic-gradient-maxLik.Rnw:575-584
###################################################
gradloss <- function(theta, index)  {
   e <- yt[index] - Xt[index,,drop=FALSE] %*% theta
   g <- -2*t(e) %*% Xt[index,,drop=FALSE]
   -g/length(index)
}
loss <- function(theta, index) {
   e <- yv - Xv %*% theta
   -crossprod(e)/length(yv)
}


###################################################
### code chunk number 14: batch1
###################################################
res <- maxSGA(loss, gradloss,
              start=start,
              nObs=nrow(Xt),  # note: only training data now
              control=list(iterlim=100,
                           SG_batchSize=1,
                           SG_learningRate=1e-5,
                           SG_clip=1e4,
                           storeParameters=TRUE,
                           storeValues=TRUE
                           )  
              )
par <- storedParameters(res)
val <- storedValues(res)
par(mfrow=c(1,2))
plot(par[,1], par[,2], type="b", pch=".",
     xlab=names(start)[1], ylab=names(start)[2], main="Parameters")
iB <- c(40, nrow(par)/2, nrow(par))
iA <- iB - 1
arrows(par[iA,1], par[iA,2], par[iB,1], par[iB,2], length=0.1)
plot(seq(length=length(val))-1, -val, type="l",
     xlab="epoch", ylab="MSE", main="Loss",
     log="y")


###################################################
### code chunk number 15: momentum
###################################################
res <- maxSGA(loss, gradloss,
              start=start,
              nObs=nrow(Xt),
              control=list(iterlim=100,
                           SG_batchSize=1,
                           SG_learningRate=1e-6,
                           SG_clip=1e4,
                           SGA_momentum = 0.99,
                           storeParameters=TRUE,
                           storeValues=TRUE
                           )  
              )
par <- storedParameters(res)
val <- storedValues(res)
par(mfrow=c(1,2))
plot(par[,1], par[,2], type="b", pch=".",
     xlab=names(start)[1], ylab=names(start)[2], main="Parameters")
iB <- c(40, nrow(par)/2, nrow(par))
iA <- iB - 1
arrows(par[iA,1], par[iA,2], par[iB,1], par[iB,2], length=0.1)
plot(seq(length=length(val))-1, -val, type="l",
     xlab="epoch", ylab="MSE", main="Loss",
     log="y")


###################################################
### code chunk number 16: Adam
###################################################
res <- maxAdam(loss, gradloss,
              start=start,
              nObs=nrow(Xt),
              control=list(iterlim=100,
                           SG_batchSize=1,
                           SG_learningRate=1e-6,
                           SG_clip=1e4,
                           storeParameters=TRUE,
                           storeValues=TRUE
                           )  
              )
par <- storedParameters(res)
val <- storedValues(res)
par(mfrow=c(1,2))
plot(par[,1], par[,2], type="b", pch=".",
     xlab=names(start)[1], ylab=names(start)[2], main="Parameters")
iB <- c(40, nrow(par)/2, nrow(par))
iA <- iB - 1
arrows(par[iA,1], par[iA,2], par[iB,1], par[iB,2], length=0.1)
plot(seq(length=length(val))-1, -val, type="l",
     xlab="epoch", ylab="MSE", main="Loss",
     log="y")


###################################################
### code chunk number 17: SANN
###################################################
val <- NULL
# loop over batch sizes
for(B in c(1,10,100)) {
   res <- maxAdam(loss, gradloss,
                  start=start,
                  nObs=nrow(Xt),
                  control=list(iterlim=200,
                               SG_batchSize=1,
                               SG_learningRate=1e-6,
                               SG_clip=1e4,
                               SG_patience=5,
                           # worse value allowed only 5 times
                               storeValues=TRUE
                               )  
                  )
   cat("Batch size", B, ",", nIter(res),
       "epochs, function value", maxValue(res), "\n")
   val <- c(val, na.omit(storedValues(res)))
   start <- coef(res)
}
plot(seq(length=length(val))-1, -val, type="l",
     xlab="epoch", ylab="MSE", main="Loss",
     log="y")
summary(res)


