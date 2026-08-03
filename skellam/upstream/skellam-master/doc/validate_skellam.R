## ----loading, echo = FALSE, results = 'hide', message=FALSE, warning=FALSE----
library(skellam)

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----variables, message=FALSE-------------------------------------------------
N <- 5000
lambda1 <- 1.5
lambda2 <- 0.5

## ----distribution, message=FALSE----------------------------------------------
X <- rpois(N, lambda1)
Y <- rpois(N, lambda2)
XminusY <- X - Y
Z <- rskellam(N, lambda1, lambda2)

## ----figure1, message=FALSE, fig.width=6, fig.height=4, out.width="100%"------
# Your plot code here
# Density comparison
plot(table(XminusY), xlab = "X - Y", ylab = "", type = "p", pch = 1)
points(table(Z), col = "red", type = "p", pch = 3, cex = 2)
xseq <- seq(floor(par("usr")[1]), ceiling(par("usr")[2]))
points(xseq, N * dskellam(xseq, lambda1, lambda2), col = "blue",
       pch = 4, cex = 3)
legend("topright", pch = c(1, 3, 4), col = c("black", "red", "blue"),
       legend = c("rpois-rpois", "rskellam", "dskellam"))

## ----qZ, message=FALSE, fig.width=6, fig.height=4, out.width="100%"-----------
# Quantile comparison
Sprob <- seq(0, 1, by = 1/100)
qZ <- quantile(Z, prob = Sprob)
plot(qZ, qskellam(Sprob, lambda1, lambda2))
abline(0, 1, col = "#FF000040")

## ----variables2, message=FALSE------------------------------------------------
lambda1 <- 12
lambda2 <- 8

X <- rpois(N, lambda1)
Y <- rpois(N, lambda2)
XminusY <- X - Y
Z <- rskellam(N, lambda1, lambda2)

## ----density2, message=FALSE, fig.width=6, fig.height=4, out.width="100%"-----
# Density comparison
plot(table(XminusY), xlab = "X - Y", ylab = "", type = "p", pch = 1)
points(table(Z), col = "red", type = "p", pch = 3, cex = 2)
xseq <- seq(floor(par("usr")[1]), ceiling(par("usr")[2]))
points(xseq, N * dskellam(xseq, lambda1, lambda2), col = "blue",
       pch = 4, cex = 3)
legend("topright", pch = c(1, 3, 4), col = c("black", "red", "blue"),
       legend = c("rpois-rpois", "rskellam", "dskellam"))

## ----qZ2, message=FALSE, fig.width=6, fig.height=4, out.width="100%"----------
# Quantile comparison
Sprob <- seq(0, 1, by = 1/100)
qZ <- quantile(Z, prob = Sprob)
plot(qZ, qskellam(Sprob, lambda1, lambda2))
abline(0, 1, col = "#FF000040")

