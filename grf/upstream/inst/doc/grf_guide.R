## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE,
  fig.width = 5,
  fig.height = 4,
  fig.align = "center"
)
set.seed(123)
options(digits = 2)

library(grf)

has.policytree = requireNamespace("policytree", quietly = TRUE)

## ----echo=FALSE, fig.align='center', out.height='1.375in', results='asis', out.extra='keepaspectratio'----
if (knitr::is_html_output()) {
  cat('<center><img src="https://raw.githubusercontent.com/grf-labs/grf/master/images/logo/grf_logo_green.png" height="132"></center>')
} else {
  knitr::include_graphics("grf_logo_green.png")
}

## ----echo=FALSE---------------------------------------------------------------
set.seed(2)
n = 20; p = 1
X = matrix(rnorm(n * p), n, p)
Y = X[, 1] + rnorm(n)
rf = regression_forest(X, Y, num.trees = 1, ci.group.size = 1, honesty = FALSE, sample.fraction = 1, num.threads = 1)
t = get_tree(rf, 1)
sval = t$nodes[[1]]$split_value

plot(c(X, 1.5), c(Y, 2),
     ylab = "Y",
     xlab = expression(X[j]),
     col = c(as.integer(X[,1] <= sval) + 1, 4),
     pch = c(rep(1, n), 3))
abline(v = sval)

## -----------------------------------------------------------------------------
data("schoolrct")
Y <- schoolrct$outcome
W <- schoolrct$treatment
school <- schoolrct$school
X <- schoolrct[-(1:3)]

## -----------------------------------------------------------------------------
colnames(X)

## -----------------------------------------------------------------------------
c.forest <- causal_forest(X, Y, W, W.hat = 0.5, clusters = school, num.threads = 2)

## -----------------------------------------------------------------------------
ate <- average_treatment_effect(c.forest)
ate
ate[1] / sd(Y) # effect in SD units

## -----------------------------------------------------------------------------
best_linear_projection(c.forest, X[c("is.female",
                                     "family.receives.cash.transfer",
                                     "is.unemployed",
                                     "financial.autonomy.index")])

## -----------------------------------------------------------------------------
tau.hat.oob <- predict(c.forest)$predictions
hist(tau.hat.oob)

## -----------------------------------------------------------------------------
rate.oob <- rank_average_treatment_effect(c.forest, tau.hat.oob)
t.stat.oob <- rate.oob$estimate / rate.oob$std.err
pnorm(t.stat.oob, lower.tail = FALSE)  # one-sided p-value Pr(>t)

## -----------------------------------------------------------------------------
var.imp <- variable_importance(c.forest)
head(colnames(X)[order(var.imp, decreasing = TRUE)])

## -----------------------------------------------------------------------------
rate.autonomy <- rank_average_treatment_effect(
  c.forest,
  -1 * X$financial.autonomy.index, # Negate priorities to rank by lowest first.
  q = seq(0.05, 1, length.out = 100),
  subset = !is.na(X$financial.autonomy.index))
plot(rate.autonomy, main = "TOC: By increasing financial autonomy")

## -----------------------------------------------------------------------------
rate.autonomy
rate.autonomy$estimate / rate.autonomy$std.err

## ----eval=has.policytree------------------------------------------------------
library(policytree)

Gamma.hat <- double_robust_scores(c.forest)
no.missing <- complete.cases(X)
X.nm <- X[no.missing, ]
Gamma.nm <- Gamma.hat[no.missing, ]

## ----eval=has.policytree------------------------------------------------------
colMeans(Gamma.nm)

## ----eval=has.policytree------------------------------------------------------
Gamma.nm[, "treated"] <- Gamma.nm[, "treated"] - ate["estimate"]

## ----eval=has.policytree------------------------------------------------------
ptree <- policy_tree(X.nm, Gamma.nm, depth = 1)
print(ptree)

## -----------------------------------------------------------------------------
data("attentionrct")
Y <- attentionrct$outcome.correct.ans.per.second
W <- 1 - attentionrct$treatment
X <- attentionrct[-(1:4)]

## -----------------------------------------------------------------------------
train <- sample(nrow(X), nrow(X) * 0.6)
test <- -train

Y.forest <- regression_forest(X[train, ], Y[train], num.trees = 500, num.threads = 2)
varimp.Y <- variable_importance(Y.forest)
keep <- colnames(X)[order(varimp.Y, decreasing = TRUE)[1:10]]
X.cate <- X[, keep]
print(keep)

## -----------------------------------------------------------------------------
# Estimate CATE function.
Y.hat.train <- predict(Y.forest)$predictions
cate.forest <- causal_forest(
  X.cate[train, ],
  Y[train],
  W[train],
  Y.hat = Y.hat.train,
  W.hat = 0.5,
  num.threads = 2
)
tau.hat.eval <- predict(cate.forest, X.cate[test, ])$predictions

# Fit evaluation forest.
eval.forest <- causal_forest(
  X[test, ],
  Y[test],
  W[test],
  W.hat = 0.5,
  num.threads = 2
)

## -----------------------------------------------------------------------------
rate.cate <- rank_average_treatment_effect(
  eval.forest,
  tau.hat.eval,
  q = seq(0.05, 1, length.out = 100)
)
print(rate.cate)
rate.cate$estimate / rate.cate$std.err # t-statistic

plot(rate.cate, main = "TOC: By decreasing CATEs")

## ----echo=FALSE, fig.align='center', out.height='0.667in', results='asis', out.extra='keepaspectratio'----
if (knitr::is_html_output()) {
  cat('<center><img src="https://raw.githubusercontent.com/grf-labs/grf/master/images/logo/grf_leaf_green.png" height="64"></center>')
} else {
  knitr::include_graphics("grf_leaf_green.png")
}

