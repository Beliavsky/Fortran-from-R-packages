#' Calculate purchasing power
#' 
#' Calculates changes in purchasing power over time,
#' taking into account the real interest rate.
#' 
#' The real interest rate is the interest rate after inflation.
#' If negative (e.g. equal to the average yearly inflation rate) 
#' it can show diminishing purchasing power over time.
#' If positive, it can show increasing purchasing power over time,
#' and effect of compounding interest on the purchasing power.
#' 
#' @seealso 
#' \itemize{
#'   \item [How to Determine Our Optimal Asset Allocation?](https://www.r4good.academy/en/blog/optimal-asset-allocation/index.en.html#why-keeping-all-your-savings-in-cash-isnt-the-best-idea)
#' }
#' 
#' @param x A numeric. The initial amount of money.
#' @param years A numeric. The number of years.
#' @param real_interest_rate A numeric. The yearly real interest rate.
#' 
#' @return A numeric. The purchasing power.
#' 
#' @examples
#' calc_purchasing_power(x = 10, years = 30, real_interest_rate = -0.02)
#' calc_purchasing_power(x = 10, years = 30, real_interest_rate = 0.02)
#' 
#' @export
#' 
calc_purchasing_power <- function(x, years, real_interest_rate) {
  
  purchasing_power <- ifelse(real_interest_rate < 0, 
                             x / (1 + abs(real_interest_rate)) ^ years, 
                             x * (1 + real_interest_rate) ^ years)
  
  return(purchasing_power)
}
