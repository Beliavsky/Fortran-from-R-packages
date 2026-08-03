#' Calculating retirement ruin probability
#' 
#' @param portfolio_return_mean A numeric. Mean of portfolio returns.
#' @param portfolio_return_sd A numeric. Standard deviation of portfolio returns.
#' @param age A numeric. Current age.
#' @param gompertz_mode A numeric. Gompertz mode.
#' @param gompertz_dispersion A numeric. Gompertz dispersion.
#' @param portfolio_value A numeric. Initial portfolio value.
#' @param monthly_spendings A numeric. Monthly spendings.
#' @param yearly_spendings A numeric. Yearly spendings.
#' @param spending_rate A numeric. Spending rate 
#' (initial withdrawal rate).
#' 
#' @return A numeric. The probability of retirement ruin (between 0 and 1),
#' representing the likelihood of running out of money during retirement.
#' 
#' @references Milevsky, M.A. (2020). Retirement Income Recipes in R: From Ruin Probabilities to Intelligent Drawdowns. Use R! Series. \doi{10.1007/978-3-030-51434-1}.
#' 
#' @examples
#' calc_retirement_ruin(
#'   age                   = 65,
#'   gompertz_mode         = 88,
#'   gompertz_dispersion   = 10,
#'   portfolio_value       = 1000000,
#'   monthly_spendings     = 3000,  
#'   portfolio_return_mean = 0.02,
#'   portfolio_return_sd   = 0.15
#' )
#' @export

calc_retirement_ruin <- function(
  portfolio_return_mean,
  portfolio_return_sd,
  age,
  gompertz_mode,
  gompertz_dispersion,
  portfolio_value,
  monthly_spendings,
  yearly_spendings = 12 * monthly_spendings,
  spending_rate    = yearly_spendings / portfolio_value
) {

  nu    <- portfolio_return_mean
  sigma <- portfolio_return_sd
  x     <- age
  m     <- gompertz_mode
  b     <- gompertz_dispersion

  mu <- nu + (0.5) * sigma^2 
  M1 <- calc_a(mu - sigma^2, x, m, b)
  M2 <- (M1 - calc_a(2 * mu - 3 * sigma^2, x, m, b)) / (mu / 2 - sigma^2)
  
  alpha <- (2 * M2 - M1^2) / (M2 - M1^2)
  beta  <- (M2 - M1^2) / (M2 * M1)

  stats::pgamma(
    q          = spending_rate,
    shape      = alpha,
    scale      = beta,
    lower.tail = TRUE
  )
}

calc_a <- function(v, x, m, b) {
  
  b * exp(exp((x - m) / b) + (x - m) * v) * 
    calc_incomplete_gamma(-b * v, exp((x - m) / b))
}

calc_incomplete_gamma <- function(a, c) {

  integrand <- function(t) {
    t^(a - 1) * exp(-t)
  }

  stats::integrate(integrand, c, Inf)$value
}
