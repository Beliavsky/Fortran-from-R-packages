#' Calculate optimal risky asset allocation
#' 
#' Calculates the optimal allocation to the risky asset 
#' using the Merton Share formula.
#' 
#' Can be used to calculate the optimal allocation to the risky asset
#' for vectors of inputs.
#' 
#' @seealso 
#' * [How to Determine Our Optimal Asset Allocation?](https://www.r4good.academy/en/blog/optimal-asset-allocation/index.en.html#what-do-you-need-to-calculate-your-optimal-asset-allocation)
#' * Haghani V., White J. (2023) "The Missing Billionaires: A Guide to Better Financial Decisions." ISBN:978-1-119-74791-8.
#' 
#' @param risky_asset_return_mean A numeric. 
#' The expected (average) yearly return of the risky asset.
#' @param risky_asset_return_sd A numeric. 
#' The standard deviation of the yearly returns of the risky asset.
#' @param safe_asset_return A numeric. 
#' The expected yearly return of the safe asset.
#' @param risk_aversion A numeric. 
#' The risk aversion coefficient.
#' 
#' @returns A numeric. 
#' The optimal allocation to the risky asset.
#' In case of [NaN()] (because of division by zero) 
#' the optimal allocation to the risky asset is set to 0.
#' 
#' @examples
#' calc_optimal_risky_asset_allocation(
#'   risky_asset_return_mean = 0.05,
#'   risky_asset_return_sd   = 0.15,
#'   safe_asset_return       = 0.02,
#'   risk_aversion           = 2
#' )
#' 
#' calc_optimal_risky_asset_allocation(
#'   risky_asset_return_mean = c(0.05, 0.06),
#'   risky_asset_return_sd   = c(0.15, 0.16),
#'   safe_asset_return       = 0.02,
#'   risk_aversion           = 2
#' )
#' @export
#' 
calc_optimal_risky_asset_allocation <- function(risky_asset_return_mean,
                                                risky_asset_return_sd,
                                                safe_asset_return,
                                                risk_aversion) {
  
  risky_asset_excess_return <- risky_asset_return_mean - safe_asset_return

  optimal_risky_asset_allocation <-
    risky_asset_excess_return / (risk_aversion * risky_asset_return_sd ^ 2)

  optimal_risky_asset_allocation[is.nan(optimal_risky_asset_allocation)] <- 0

  optimal_risky_asset_allocation
}
