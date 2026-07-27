#'Computes the empirical Conditional Value-at-Risk, Value-at-Risk and Mean Absolute Deviation for losses with given probabilities
#'
#'@param returns vector of losses
#'@param prob vector of probability of losses
#'@param alpha confidence level for CVaR and VaR
#'@return list of values
#'@keywords internal 

.RISK_post <- function (returns, prob, alpha)
{
  rdigits = 10
  o = order (returns)
  sorted_returns = as.matrix(returns[o])
  weight = as.matrix(prob[o])
  mu = (t(sorted_returns) %*% weight)[1,1]
  mu = matrix(rep(mu,length(o)), length(o),1)
  index = sum(cumsum(weight) < alpha) +1

  tail_sum =  t(sorted_returns[(index + 1):length(o)]) %*% (weight[(index + 1):length(o)])  + sorted_returns [index] * (sum(weight[1:index]) - alpha)
  CVaR = as.numeric(tail_sum / (1 - alpha))
  MAD =  abs(t(sorted_returns -mu)) %*% weight

  return (list (var = round(sorted_returns [index], digits = rdigits), cvar = round(CVaR, digits = rdigits), mad = round(MAD, digits = rdigits)))
}
