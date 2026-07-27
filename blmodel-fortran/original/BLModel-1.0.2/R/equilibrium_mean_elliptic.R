#'Solves the inverse portfolio optimization problem with risk measured by variance to find equilibrium returns
#'
#'@param MCov g
#'@param MarketPortfolio g
#'@param MarketPriceOfRisk g
#'@return res list contains market_returns and portfolio
#'
#'

.equilibrium_mean_elliptic  <- function (MCov, MarketPortfolio, MarketPriceOfRisk)
{
  Portfolio = cbind(MarketPortfolio)
  if( is.null(rownames(MCov))) {
    clab <- as.character(1:nrow(MCov))
  } else {
    clab <- rownames(MCov)
  }

  rownames(Portfolio) = clab
  colnames(Portfolio)= as.character("assets weigths")

  mu = MarketPriceOfRisk * MCov %*%  Portfolio
  rownames(mu) = clab
  colnames(mu)= as.character("assets excess returns")
  res = list(market_returns = mu, portfolio = Portfolio)
  return (res)
}
