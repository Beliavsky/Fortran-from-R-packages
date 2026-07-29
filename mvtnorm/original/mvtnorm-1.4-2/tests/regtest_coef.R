
library("mvtnorm")

options(digits = 4)

J <- 5

set.seed(290875)

lxd <- ltMatrices(runif(J * (J + 1) / 2) + 1, diag = TRUE, names = LETTERS[seq_len(J)])
m <- matrix(1:J, ncol = 1)
d <- Mult(lxd, m)

d0 <- mvnorm(invchol = lxd)
Y <- as.data.frame(t(simulate(d0, nsim = 1e6, seed = 29)))
mlm <- lm(D ~ 0 + A + B + C, data = Y)
(cf <- coef(d0, which = "D"))
coef(mlm)
summary(mlm)$sigma
md <- margDist(d0, which = 1:4)
cd <- condDist(md, which_given = 1:3, given = diag(3))
class(cd$scale)
coef(cd)

d0 <- mvnorm(chol = solve(lxd))
Y <- as.data.frame(t(simulate(d0, nsim = 1e6, seed = 29)))
mlm <- lm(D ~ 0 + A + B + C, data = Y)
(cf <- coef(d0, which = "D"))
coef(mlm)
summary(mlm)$sigma
md <- margDist(d0, which = 1:4)
cd <- condDist(md, which_given = 1:3, given = diag(3))
coef(cd)

(d1 <- mvnorm(invchol = lxd, mean = m))
Y <- as.data.frame(t(simulate(d1, nsim = 1e6, seed = 29)))
mlm <- lm(D ~ A + B + C, data = Y)
(cf <- coef(d1, which = "D"))
coef(mlm)
summary(mlm)$sigma
md <- margDist(d1, which = 1:4)

cd <- condDist(md, which_given = 1:3, given = diag(3))
class(cd$scale)
coef(cd)
cf[1] + cf[-1]

d2 <- mvnorm(invchol = lxd, invcholmean = d)
Y <- as.data.frame(t(simulate(d2, nsim = 1e6, seed = 29)))
mlm <- lm(D ~ A + B + C, data = Y)
coef(d2, which = "D")
coef(mlm)
summary(mlm)$sigma
coef(condDist(margDist(d2, which = 1:4), which_given = 1:3, given = diag(3)))
cf[1] + cf[-1]

d3 <- mvnorm(chol = solve(lxd), invcholmean = d)
Y <- as.data.frame(t(simulate(d3, nsim = 1e6, seed = 29)))
mlm <- lm(D ~ A + B + C, data = Y)
coef(d3, which = "D")
coef(mlm)
summary(mlm)$sigma
coef(condDist(margDist(d3, which = 1:4), which_given = 1:3, given = diag(3)))
cf[1] + cf[-1]

d4 <- mvnorm(chol = solve(lxd), mean = m)
Y <- as.data.frame(t(simulate(d4, nsim = 1e6, seed = 29)))
mlm <- lm(D ~ A + B + C, data = Y)
coef(d4, which = "D")
coef(mlm)
summary(mlm)$sigma
coef(condDist(margDist(d4, which = 1:4), which_given = 1:3, given = diag(3)))
cf[1] + cf[-1]

