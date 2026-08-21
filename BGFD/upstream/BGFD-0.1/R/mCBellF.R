#' @export
#' @import AdequacyModel
#' @import graphics
mCBellF<-function (x, a, b, lambda, method="B")
{

pdf_CBellF<-function(par,x){
a=par[1]
b=par[2]
lambda=par[3]
G=(1-((1+(x/a)^b))^(-1))
g=((1+(x/a)^b))^(-2) *b*x^(b-1)*a^(-b)
F=(exp(exp(lambda*G)-1)-1)/(exp(exp(lambda)-1)-1)
f=lambda*g*exp(lambda*G)*exp(exp(lambda*G)-1)/(exp(exp(lambda)-1)-1)
return(f)
}
cdf_CBellF<-function(par,x){
a=par[1]
b=par[2]
lambda=par[3]
  G=(1-((1+(x/a)^b))^(-1))
g=((1+(x/a)^b))^(-2) *b*x^(b-1)*a^(-b)
F=(exp(exp(lambda*G)-1)-1)/(exp(exp(lambda)-1)-1)
f=lambda*g*exp(lambda*G)*exp(exp(lambda*G)-1)/(exp(exp(lambda)-1)-1)
return(F)
}
res = suppressWarnings(AdequacyModel::goodness.fit(pdf = pdf_CBellF, cdf = cdf_CBellF, starts = c(a, b, lambda), data = x, method = method, mle = NULL))
aux = cbind(res$mle, res$Erro)
colnames(aux) = c("MLE", "SE")
aux1 = cbind(res$AIC, res$BIC, res$W,res$A, res$Value)
  colnames(aux1) = c("AIC",  "BIC",  "W", "A","-2L")
  rownames(aux1) = c("")
  aux2 = cbind(res$KS$statistic, res$KS$p.value)
  colnames(aux2) = c("KS Statistic", "KS p-value")
  rownames(aux2) = c("")
  aux3 = cbind(if (res$Convergence == 0) {
    "Converged"
  }
  else {
    "Not Converged"
  })

	  
  colnames(aux3) = c("")
  rownames(aux3) = c("")
  list(Estimates = aux, `Goodness-of-Fit Tests`  = aux1, `Kolmogorov-Smirnov Test` = aux2,
       `Convergence Status` = aux3)
  }
