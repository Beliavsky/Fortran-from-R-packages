durbinH <- function (model,Ly.label="ar1") {
  d <- car::durbinWatsonTest(as.vector(resid(model)))
  n <- length(fitted(model)) -length(coef(model))+ 1
  v <- diag(vcov(model))[Ly.label]
  durbinH <- as.numeric((1 - d/2) * sqrt(n/(1 - n * v)))
  return(list(Durbin.h=durbinH,pvalue=(1-pnorm(abs(durbinH)))*2))
}