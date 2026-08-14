
library("mvtnorm")

options(warn = -1L)

chk <- function(...) isTRUE(all.equal(...))

set.seed(29)

J <- 500
N <- 10

X <- diag(J)
X[lower.tri(X, diag = TRUE)] <- prm <- runif(J * (J + 1) / 2)
L <- ltMatrices(prm, diag = TRUE, byrow = FALSE, names = FALSE)

Y <- matrix(rnorm(J * N), nrow = J)

object.size(X)
object.size(prm)
object.size(L)

system.time(a <- X %*% Y)
system.time(b <- L %*% Y)

chk(a, b, check.attributes = FALSE)

system.time(a <- crossprod(X, Y))
system.time(b <- crossprod(L, Y))

chk(a, b, check.attributes = FALSE)

system.time(a <- crossprod(X))
system.time(b <- crossprod(L))

chk(a, as.array(b)[,,1], check.attributes = FALSE)

system.time(a <- tcrossprod(X))
system.time(b <- tcrossprod(L))
system.time(b1 <- vectrick(L, transpose = c(FALSE, TRUE)))

chk(a, as.array(b)[,,1], check.attributes = FALSE)
chk(unclass(b), unclass(b1), check.attributes = FALSE)

system.time(L <- ltMatrices(L, byrow = TRUE))

system.time(a <- X %*% Y)
system.time(b <- L %*% Y)

chk(a, b, check.attributes = FALSE)

system.time(a <- crossprod(X, Y))
system.time(b <- crossprod(L, Y))

chk(a, b, check.attributes = FALSE)

system.time(a <- crossprod(X))
system.time(b <- crossprod(L))

chk(a, as.array(b)[,,1], check.attributes = FALSE)

system.time(a <- tcrossprod(X))
system.time(b <- tcrossprod(L))

chk(a, as.array(b)[,,1], check.attributes = FALSE)

