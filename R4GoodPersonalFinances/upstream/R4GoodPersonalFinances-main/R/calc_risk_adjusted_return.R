#' Calculate risk adjusted return
#' 
#' Calculates the risk adjusted return for portfolio
#' of given allocation to the risky asset.
#' 
#' @inheritParams calc_optimal_risky_asset_allocation
#' 
#' @seealso 
#' * [How to Determine Our Optimal Asset Allocation?](https://www.r4good.academy/en/blog/optimal-asset-allocation/index.en.html#how-much-risk-is-enough)
#' * Haghani V., White J. (2023) "The Missing Billionaires: A Guide to Better Financial Decisions." ISBN:978-1-119-74791-8.
#' 
#' @param risky_asset_allocation A numeric. 
#' The allocation to the risky asset. Could be a vector.
#' If it is the optimal allocation then parameters
#' `risky_asset_return_sd` and `risk_aversion`
#' can be omitted.
#' 
#' @returns A numeric. The risk adjusted return.
#' 
#' @examples
#' calc_risk_adjusted_return(
#'   safe_asset_return = 0.02,
#'   risky_asset_return_mean = 0.04,
#'   risky_asset_return_sd = 0.15,
#'   risky_asset_allocation = 0.5,
#'   risk_aversion = 2
#' )
#' 
#' calc_risk_adjusted_return(
#'   safe_asset_return = 0.02,
#'   risky_asset_return_mean = 0.04,
#'   risky_asset_allocation = c(0.25, 0.5, 0.75),
#'   risky_asset_return_sd = 0.15,
#'   risk_aversion = 2
#' )
#' @export
#' 
calc_risk_adjusted_return <- function(safe_asset_return,
                                      risky_asset_return_mean,
                                      risky_asset_allocation,
                                      risky_asset_return_sd = NULL,risk_aversion         = NULL
                                    ) {
  
  if (
    (!is.null(risky_asset_return_sd) && is.null(risk_aversion)) ||
    (is.null(risky_asset_return_sd) && !is.null(risk_aversion))
    ) {
      stop(
        "Both `risky_asset_return_sd` and `risk_aversion` cannot be NULL.",
        call. = FALSE
      )
    }
  
  risky_asset_excess_return <- risky_asset_return_mean - safe_asset_return
                                      
  if (is.null(risky_asset_return_sd) & is.null(risk_aversion)) {

    return(
      safe_asset_return +
        (risky_asset_allocation * risky_asset_excess_return) / 2
    )
  }

  safe_asset_return + risky_asset_allocation * (
    risky_asset_excess_return - (
      (risky_asset_allocation * risk_aversion * risky_asset_return_sd ^ 2) / 2
    )
  )
  }
