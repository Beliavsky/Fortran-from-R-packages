## -----------------------------------------------------------------------------
#| echo: false
#| warning: false
library(COMPoissonReg)
set.seed(1235)


## -----------------------------------------------------------------------------
control = get.control(
	ymax = 100000,
	hybrid.tol = 1e-2,
	truncate.tol = 1e-6
)


## -----------------------------------------------------------------------------
control = getOption("COMPoissonReg.control")
control$ymax
control$hybrid.tol
control$truncate.tol


## -----------------------------------------------------------------------------
options(COMPoissonReg.control = control)


## -----------------------------------------------------------------------------
ncmp(lambda = 1.5, nu = 1.2)
ncmp(lambda = 1.5, nu = 1.2, log = TRUE)
ncmp(lambda = 1.5, nu = 1.2, log = TRUE, control = get.control(hybrid.tol = 1e10))
ncmp(lambda = 1.5, nu = 1.2, log = TRUE, control = get.control(hybrid.tol = 1e-10))


## -----------------------------------------------------------------------------
print_warning = function(x) { print(strwrap(x), quote = FALSE) }


## -----------------------------------------------------------------------------
nu_seq = c(1, 0.5, 0.2, 0.1, 0.05, 0.03)
tryCatch({ tcmp(lambda = 1.5, nu = nu_seq) }, warning = print_warning)


## -----------------------------------------------------------------------------
tcmp(lambda = 1.5, nu = nu_seq, control = get.control(ymax = 3e6))


## -----------------------------------------------------------------------------
tcmp(lambda = 1.2, nu = 0.03, control = get.control(ymax = 1200))


## -----------------------------------------------------------------------------
#| prompt: false
#| warning: false
library(ggplot2)

nu_seq = seq(0.03, 1.5, length.out = 20)
nc1 = ncmp(lambda = 0.5, nu = nu_seq, log = TRUE)
nc2 = ncmp(lambda = 1.05, nu = nu_seq, log = TRUE)
nc3 = ncmp(lambda = 1.20, nu = nu_seq, log = TRUE)


## -----------------------------------------------------------------------------
#| fig-width: 5
#| fig-height: 3
#| fig-align: "center"
#| prompt: false
#| fig-cap: "Log of normalizing constant for $\\lambda = 0.5$ ($\\circ$), $\\lambda = 1.05$ ($\\Delta$), and $\\lambda = 1.20$ ($+$)."

ggplot() +
	geom_point(data = data.frame(x = nu_seq, y = nc1), aes(x = x, y = y), pch = 1) +
	geom_point(data = data.frame(x = nu_seq, y = nc2), aes(x = x, y = y), pch = 2) +
	geom_point(data = data.frame(x = nu_seq, y = nc3), aes(x = x, y = y), pch = 3) +
	xlab("nu") +
	ylab("log of normalizing constant") +
	theme_bw()


## -----------------------------------------------------------------------------
dcmp(0, lambda = 10, nu = 0.9)
dcmp(0:17, lambda = 10, nu = 0.9, log = TRUE)
dcmp(c(0, 1, 2), lambda = c(10, 11, 12), nu = c(0.9, 1.0, 1.1), log = TRUE)


## -----------------------------------------------------------------------------
rcmp(50, lambda = 10, nu = 0.9)


## -----------------------------------------------------------------------------
pcmp(0:17, lambda = 10, nu = 0.9)


## -----------------------------------------------------------------------------
qq = seq(0, 0.95, length.out = 10)
qcmp(qq, lambda = 10, nu = 0.9)


## -----------------------------------------------------------------------------
tryCatch({ rcmp(1, lambda = 2, nu = 0.01) }, warning = print_warning)


## -----------------------------------------------------------------------------
tryCatch({
	qcmp(0.9999999, lambda = 1.5, nu = 0.5)
}, warning = print_warning)


## -----------------------------------------------------------------------------
#| fig-width: 3
#| fig-height: 3
#| fig-align: "center"
#| prompt: false
#| fig-show: "hold"

library(ggplot2)

n = 100000
lambda = 0.5
nu = 0.1
x = rcmp(n, lambda, nu)

xx = seq(-1, max(x))  ## Include -1 to ensure it gets probability zero
qq = seq(0, 0.99, length.out = 100)

fx = dcmp(xx, lambda, nu)
px = pcmp(xx, lambda, nu)
qx = qcmp(qq, lambda, nu)

qx_emp = quantile(x, probs = qq)


## -----------------------------------------------------------------------------
#| fig-width: 3
#| fig-height: 3
#| prompt: false
#| fig-cap: "Empirical density of draws (histogram) versus density computed via the dcmp function (points)."
ggplot() +
	geom_bar(data = data.frame(x = x), aes(x = x, y = ..prop..), fill = "NA",
		col = "black") +
	geom_point(data = data.frame(x = xx[-1], fx = fx[-1]), aes(x, fx)) +
	ylab("Density") +
	theme_bw()


## -----------------------------------------------------------------------------
#| fig-width: 3
#| fig-height: 3
#| prompt: false
#| fig-cap: "Empirical CDF of draws (solid line) versus CDF computed via the pcmp function (points)."
ggplot() +
	stat_ecdf(data = data.frame(x = x), aes(x), geom = "step") +
	geom_point(data = data.frame(x = xx, px = px), aes(x, px)) +
	ylab("Probability") +
	theme_bw()


## -----------------------------------------------------------------------------
#| fig-width: 3
#| fig-height: 3
#| prompt: false
#| fig-cap: "Empirical quantiles of draws (`o`) versus quantiles computed via the qcmp function (`+`)."
ggplot() +
	geom_point(data = data.frame(x = qq, qx_emp = qx_emp), aes(qq, qx_emp), pch = 1) +
	geom_point(data = data.frame(x = qq, qx = qx), aes(qq, qx), pch = 3) +
	xlab("Probability") +
	ylab("Quantile") +
	theme_bw()


## -----------------------------------------------------------------------------
ecmp(lambda = 10, nu = 1.2)
ecmp(lambda = 1.5, nu = 0.5)
ecmp(lambda = 1.5, nu = 0.05)
ecmp(lambda = 1.5, nu = 0.05, control = get.control(hybrid.tol = 1e-10))
ecmp(lambda = 1.5, nu = 0.05, control = get.control(hybrid.tol = 1e10))


## -----------------------------------------------------------------------------
vcmp(lambda = 10, nu = 1.2)
vcmp(lambda = 1.5, nu = 0.5)
vcmp(lambda = 1.5, nu = 0.05)
vcmp(lambda = 1.5, nu = 0.05, control = get.control(hybrid.tol = 1e-10))
vcmp(lambda = 1.5, nu = 0.05, control = get.control(hybrid.tol = 1e10))


## -----------------------------------------------------------------------------
M = tcmp(lambda = 1.5, nu = 0.05)
print(M)
xx = seq(0, M)
sum(xx^3 * dcmp(xx, lambda, nu))    # E(X^3)
sum(xx^4 * dcmp(xx, lambda, nu))    # E(X^4)


## -----------------------------------------------------------------------------
qq = seq(0, 0.95, length.out = 20)
rzicmp(20, lambda = 1.5, nu = 0.2, p = 0.25)
dzicmp(c(0, 1, 2), lambda = 1.5, nu = 0.2, p = 0.25)
pzicmp(c(0, 1, 2), lambda = 1.5, nu = 0.2, p = 0.25)
qzicmp(qq, lambda = 1.5, nu = 0.2, p = 0.25)


## -----------------------------------------------------------------------------
tryCatch({
	qzicmp(0.9999999, lambda = 1.5, nu = 0.5, p = 0.5)
}, warning = print_warning)


## -----------------------------------------------------------------------------
#| fig-width: 3
#| fig-height: 3
#| fig-align: "center"
#| prompt: false
#| fig-show: "hold"

library(ggplot2)

n = 100000
lambda = 0.5
nu = 0.1
p = 0.5
x = rzicmp(n, lambda, nu, p)

xx = seq(-1, max(x))  ## Include -1 to ensure it gets probability zero
qq = seq(0, 0.99, length.out = 100)

fx = dzicmp(xx, lambda, nu, p)
px = pzicmp(xx, lambda, nu, p)
qx = qzicmp(qq, lambda, nu, p)

qx_emp = quantile(x, probs = qq)


## -----------------------------------------------------------------------------
#| fig-width: 3
#| fig-height: 3
#| prompt: false
#| fig-cap: "Empirical density of draws (histogram) versus density computed via the dzicmp function (points)."

ggplot() +
	geom_bar(data = data.frame(x = x), aes(x = x, y = ..prop..), fill = "NA",
		col = "black") +
	geom_point(data = data.frame(x = xx[-1], fx = fx[-1]), aes(x, fx)) +
	ylab("Density") +
	theme_bw()


## -----------------------------------------------------------------------------
#| fig-width: 3
#| fig-height: 3
#| prompt: false
#| fig-cap: "Empirical CDF of draws (solid line) versus CDF computed via the pzicmp function (points)."
ggplot() +
	stat_ecdf(data = data.frame(x = x), aes(x), geom = "step") +
	geom_point(data = data.frame(x = xx, px = px), aes(x, px)) +
	ylab("Probability") +
	theme_bw()


## -----------------------------------------------------------------------------
#| fig-width: 3
#| fig-height: 3
#| prompt: false
#| fig-cap: "Empirical quantiles of draws (`o`) versus quantiles computed via the qzicmp function (`+`)."
ggplot() +
	geom_point(data = data.frame(x = qq, qx_emp = qx_emp), aes(qq, qx_emp), pch = 1) +
	geom_point(data = data.frame(x = qq, qx = qx), aes(qq, qx), pch = 3) +
	xlab("Probability") +
	ylab("Quantile") +
	theme_bw()


## -----------------------------------------------------------------------------
ezicmp(lambda = 1.5, nu = 0.5, p = 0.1)
ezicmp(lambda = 1.5, nu = 0.5, p = c(0.1, 0.2, 0.5))


## -----------------------------------------------------------------------------
vzicmp(lambda = 1.5, nu = 0.5, p = 0.1)
vzicmp(lambda = 1.5, nu = 0.5, p = c(0.1, 0.2, 0.5))


## -----------------------------------------------------------------------------
#| eval: false
#| prompt: false
# out = glm.cmp(formula.lambda, formula.nu = ~ 1, formula.p = NULL,
# 	data = NULL, init = NULL, fixed = NULL, control = NULL, ...)


## -----------------------------------------------------------------------------
get.init(beta = c(1, 2, 3), gamma = c(-1, 1), zeta = c(-2, -1))
get.init.zero(d1 = 3, d2 = 2, d3 = 2)


## -----------------------------------------------------------------------------
get.fixed(beta = c(1L, 2L), gamma = c(1L))


## -----------------------------------------------------------------------------
#| eval: false
# model.matrix(formula.lambda, data = data)
# model.matrix(formula.nu, data = data)
# model.matrix(formula.p, data = data)


## -----------------------------------------------------------------------------
control = getOption("COMPoissonReg.control")
control$optim.method
control$optim.control


## -----------------------------------------------------------------------------
data(freight)
print(freight)


## -----------------------------------------------------------------------------
glm.out = glm(broken ~ transfers, data = freight, family = poisson)
summary(glm.out)


## -----------------------------------------------------------------------------
cmp.out = glm.cmp(broken ~ transfers, data = freight)
print(cmp.out)


## -----------------------------------------------------------------------------
cmp2.out = glm.cmp(broken ~ transfers, formula.nu = ~ transfers, data = freight)
print(cmp2.out)


## -----------------------------------------------------------------------------
control = get.control(optim.control = list(maxit = 5, trace = 3, REPORT = 1))
cmp3.out = glm.cmp(broken ~ transfers, data = freight, control = control)


## -----------------------------------------------------------------------------
#| results: 'hide'
y = freight$broken
x = freight$transfers
glm.cmp(y ~ x)


## -----------------------------------------------------------------------------
#| results: 'hide'
freight$offx = 13
freight$offs = 1
glm.cmp(broken ~ transfers + offset(offx), data = freight)
glm.cmp(broken ~ transfers + offset(offx), formula.nu = ~1 + offset(offs), data = freight)


## -----------------------------------------------------------------------------
#| results: 'hide'
y = freight$broken
X = model.matrix(~ transfers, data = freight)
S = model.matrix(~ 1, data = freight)
offs = get.offset(x = rep(13, nrow(freight)), s = rep(1, nrow(freight)))
cmp.raw.out = glm.cmp.raw(y, X, S, offset = offs)


## -----------------------------------------------------------------------------
logLik(cmp.out)               ## Log-likelihood evaluated at MLE.
AIC(cmp.out)                  ## AIC evaluated at MLE.
BIC(cmp.out)                  ## BIC evaluated at MLE.
coef(cmp.out)                 ## Estimates of theta as a flat vector
coef(cmp.out, type = "list")  ## Estimates of theta as a named list
vcov(cmp.out)                 ## Estimated covariance matrix of theta hat
sdev(cmp.out)                 ## Standard deviations from vcov(...) diagonals
sdev(cmp.out, type = "list")  ## Standard deviations as a named list


## -----------------------------------------------------------------------------
predict(cmp.out)
predict(cmp.out, type = "link")


## -----------------------------------------------------------------------------
#| fig-width: 3
#| fig-height: 3
#| fig-align: "center"
#| prompt: false

# Prepare new data to fit by formula interface
new.df = data.frame(transfers = 0:10)

# Prepare new data to fit by raw interface
X = model.matrix(~ transfers, data = new.df)
S = model.matrix(~ 1, data = new.df)
new.data = get.modelmatrix(X = X, S = S)

# Pass new data to model from by formula interface
y.hat.new = predict(cmp.out, newdata = new.df)

# Pass new data to model from by raw interface
y.hat.new = predict(cmp.raw.out, newdata = new.data)

# Compute predictions for links
predict.out = predict(cmp.out, newdata = new.df, type = "link")

# Plot predictions
ggplot() +
	geom_point(data = new.df, aes(transfers, y.hat.new)) +
	xlab("Number of transfers") +
	ylab("Predicted number broken") +
	theme_bw()


## -----------------------------------------------------------------------------
print(y.hat.new)
print(predict.out)


## -----------------------------------------------------------------------------
leverage(cmp.out)


## -----------------------------------------------------------------------------
res.raw = residuals(cmp.out)
res.qtl = residuals(cmp.out, type = "quantile")


## -----------------------------------------------------------------------------
link.hat = predict(cmp.out, type = "link")
vv = vcmp(link.hat$lambda, link.hat$nu)
hh = leverage(cmp.out)
res.pearson = res.raw / sqrt(vv*(1-hh))


## -----------------------------------------------------------------------------
#| fig-width: 3
#| fig-height: 3
#| fig-align: "center"
#| prompt: false
#| fig-show: "hold"
#| layout-ncol: 2

plot.fit.res = function(y.hat, res) {
	ggplot(data.frame(y = y.hat, res = res)) +
		geom_point(aes(y, res)) +
		xlab("Fitted Value") +
		ylab("Residual Value") +
		theme_bw() +
		theme(plot.title = element_text(size = 10))
}

plot.qq.res = function(res) {
	ggplot(data.frame(res = res), aes(sample = res)) +
		stat_qq() +
		stat_qq_line() +
		theme_bw() +
		theme(plot.title = element_text(size = 10))
}

y.hat = predict(cmp.out)

plot.fit.res(y.hat, res.raw) +
	ggtitle("Fitted Values vs. Raw Residuals")
plot.qq.res(res.raw) +
	ggtitle("Q-Q Plot of Raw Residuals")

plot.fit.res(y.hat, res.pearson) +
	ggtitle("Fitted Values vs. Pearson Residuals")
plot.qq.res(res.pearson) +
	ggtitle("Q-Q Plot of Pearson Residuals")

plot.fit.res(y.hat, res.qtl) +
	ggtitle("Fitted Values vs. Quantile Residuals")
plot.qq.res(res.qtl) +
	ggtitle("Q-Q Plot of Quantile Residuals")


## -----------------------------------------------------------------------------
mean(res.raw^2)


## -----------------------------------------------------------------------------
equitest(cmp.out)


## -----------------------------------------------------------------------------
deviance(cmp.out)


## -----------------------------------------------------------------------------
cmp.boot = parametric.bootstrap(cmp.out, reps = 100)
head(cmp.boot)


## -----------------------------------------------------------------------------
out = apply(cmp.boot, 2, quantile, probs = c(0.025,0.975), na.rm = TRUE)
t(out)


## -----------------------------------------------------------------------------
#| prompt: false
set.seed(1234)
n = 200
x = rnorm(n, 500, 10)
X = cbind(intercept = 1, slope = x)
S = matrix(1, n, 1)
beta_true = c(-0.05, 0.05)
gamma_true = 2
lambda_true = exp(X %*% beta_true)
nu_true = exp(S %*% gamma_true)
y = rcmp(n, lambda_true, nu_true)


## -----------------------------------------------------------------------------
summary(x)
summary(y)


## -----------------------------------------------------------------------------
tryCatch({
	glm.cmp(y ~ x, formula.nu = ~ 1)
}, error = print_warning)


## -----------------------------------------------------------------------------
glm.cmp(y ~ scale(x), formula.nu = ~ 1)


## -----------------------------------------------------------------------------
glm.cmp(y ~ log(x), formula.nu = ~ 1)


## -----------------------------------------------------------------------------
control = get.control(optim.method = "BFGS", optim.control = list(maxit = 200))
suppressWarnings({
	cmp.out = glm.cmp(y ~ x, formula.nu = ~ 1, control = control)
	print(cmp.out)
})


## -----------------------------------------------------------------------------
#| prompt: false
set.seed(1234)
n = 200
x = runif(n, 1, 2)
X = cbind(intercept = 1, slope = x)
S = matrix(1, n, 1)
beta_true = c(1, 1)
gamma_true = -0.95
lambda_true = exp(X %*% beta_true)
nu_true = exp(S %*% gamma_true)
y = rcmp(n, lambda_true, nu_true)


## -----------------------------------------------------------------------------
summary(x)
summary(y)


## -----------------------------------------------------------------------------
tryCatch({
	glm.cmp(y ~ x, formula.nu = ~ 1)
}, error = print_warning)


## -----------------------------------------------------------------------------
init = get.init(beta = beta_true, gamma = gamma_true)
glm.cmp(y ~ x, formula.nu = ~ 1, init = init)


## -----------------------------------------------------------------------------
control = get.control(optim.method = "Nelder-Mead", optim.control = list(maxit = 1000))
glm.cmp(y ~ x, formula.nu = ~ 1, control = control)


## -----------------------------------------------------------------------------
data(couple)
head(couple)


## -----------------------------------------------------------------------------
glm.out = glm(UPB ~ EDUCATION + ANXIETY, data = couple, family = poisson)
summary(glm.out)


## -----------------------------------------------------------------------------
zicmp0.out = glm.cmp(UPB ~ EDUCATION + ANXIETY,
	formula.nu = ~ 1,
	formula.p = ~ EDUCATION + ANXIETY,
 	data = couple)
print(zicmp0.out)


## -----------------------------------------------------------------------------
pred.out = predict(zicmp0.out, type = "link")
summary(pred.out$lambda)


## -----------------------------------------------------------------------------
#| prompt: false
library(numDeriv)
g = function(gamma0) {
	-ncmp(lambda = exp(-0.25), nu = exp(gamma0), log = TRUE)
}
dat = data.frame(gamma0 = seq(0, -13), g = NA, d_g = NA, d2_g = NA)
for (j in 1:nrow(dat)) {
	gamma0 = dat$gamma0[j]
	dat$g[j] = g(gamma0)
	dat$d_g[j] = numDeriv::grad(func = g, x = gamma0)
	dat$d2_g[j] = numDeriv::hessian(func = g, x = gamma0)
}


## -----------------------------------------------------------------------------
print(dat)


## -----------------------------------------------------------------------------
#| prompt: false
init = coef(zicmp0.out, type = "list")
y = couple$UPB
X = model.matrix(~ EDUCATION + ANXIETY, data = couple)
S = model.matrix(~ 1, data = couple)
W = model.matrix(~ EDUCATION + ANXIETY, data = couple)
control = get.control(optim.method = "BFGS")
zicmp.out = glm.zicmp.raw(y, X, S, W,
	init = get.init(beta = c(-1,0,0), gamma = -Inf, zeta = c(-1,0,0)),
	fixed = get.fixed(gamma = 1L), control = control)


## -----------------------------------------------------------------------------
print(zicmp.out)


## -----------------------------------------------------------------------------
logLik(zicmp.out)               ## Log-likelihood evaluated at MLE
AIC(zicmp.out)                  ## AIC evaluated at MLE
BIC(zicmp.out)                  ## BIC evaluated at MLE
coef(zicmp.out)                 ## Estimates of theta as a flat vector
coef(zicmp.out, type = "list")  ## Estimates of theta as a named list
vcov(zicmp.out)                 ## Estimated covariance matrix of theta hat
sdev(zicmp.out)                 ## Standard deviations from vcov(...) diagonals
sdev(zicmp.out, type = "list")  ## Standard deviations as a named list
equitest(zicmp0.out)            ## Likelihood ratio test for H_0: gamma = 0
tryCatch({                      ## An error is thrown for model with fixed gamma
	equitest(zicmp.out)
}, error = print_warning)


## -----------------------------------------------------------------------------
y.hat = predict(zicmp.out)  ## Fitted values based on ecmp
link.hat = predict(zicmp.out, type = "link")
head(y.hat)
head(link.hat)


## -----------------------------------------------------------------------------
#| fig-width: 3
#| fig-height: 3
#| fig-align: "center"
#| prompt: false
#| fig-show: "hold"
#| layout-ncol: 2

res.raw = residuals(zicmp.out, type = "raw")
res.qtl = residuals(zicmp.out, type = "quantile")

plot.fit.res(y.hat, res.raw) +
	ggtitle("Fitted Values vs. Raw Residuals")
plot.qq.res(res.raw) +
	ggtitle("Q-Q Plot of Raw Residuals")

plot.fit.res(y.hat, res.qtl) +
	ggtitle("Fitted Values vs. Quantile Residuals")
plot.qq.res(res.qtl) +
	ggtitle("Q-Q Plot of Quantile Residuals")


## -----------------------------------------------------------------------------
#| prompt: false
new.df = data.frame(EDUCATION = round(1:20 / 20), ANXIETY = seq(-3,3, length.out = 20))
X.new = model.matrix(~ EDUCATION + ANXIETY, data = new.df)
S.new = model.matrix(~ 1, data = new.df)
W.new = model.matrix(~ EDUCATION + ANXIETY, data = new.df)
new.data = get.modelmatrix(X.new, S.new, W.new)

# For model fit using raw interface, use get.modelmatrix to prepare new design
# matrices, offsets, etc
y.hat.new = predict(zicmp.out, newdata = new.data)

# For models fit with the formula interface, pass a data.frame with the same
# structure as used in the fit.
y.hat.new = predict(zicmp0.out, newdata = new.df)


## -----------------------------------------------------------------------------
print(y.hat.new)


## -----------------------------------------------------------------------------
#| eval: false
# zicmp.boot = parametric.bootstrap(zicmp.out, reps = 100)
# head(zicmp.boot)
# apply(zicmp.boot, 2, quantile, probs = c(0.025,0.975), na.rm = TRUE)

