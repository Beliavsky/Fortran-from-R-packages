#' vamc: A package for pricing a pool of variable annuities.
#'
#' The vamc package provides a Monte Carlo engine for valuating a pool of
#' variable annuities. The key steps are:
#' YieldCurveGeneration, ScenarioGeneration, PolicyGenerationl,
#' and MonteCarloValuation.
#'
#' @section YieldCurveGeneration functions:
#' YieldCurveGeneration generates a forward curve from swap rates.
#' The forward curve is obtained by solving for swap rates that equates values of
#' floating and fixed notes.
#' @section ScenarioGeneration functions:
#' ScenarioGeneration generates a random fund scenario under Black-Scholes.
#' After simulating random index scenarios, a fundMap is used to allocate returns
#' of indices to each fund according to proportion of investment.
#' @section PolicyGenerationl functions:
#' PolicyGenerationl randomly generates a pool of variable annuities for
#' user-input birthday range, issue-date range, maturity range, account value range,
#' female percentage, fund management fee, fund base fee, product types,
#' rider fee of each type, roll-up-rate for roll-up featured guarantees,
#' withdrawal rate for GMWB, and number of policies to be generated for each type.
#' @section MonteCarloValuation functions:
#' MonteCarloValuation discounts cash flow from living and death benefits, as
#' well as risk charges for each policy in the portfolio.
#'
#' @docType package
#' @name vamc

#' @importFrom Rdpack reprompt
#' @references
#' \insertRef{gan2017vabd}{vamc}
NULL
