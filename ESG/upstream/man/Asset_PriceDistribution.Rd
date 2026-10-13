\name{Asset_PriceDistribution}
\alias{Asset_PriceDistribution}
\alias{Asset_PriceDistribution,Scenarios-method}
\title{Asset_PriceDistribution method}
\arguments{
  \item{type}{The name of the asset to price. It must be
  'Zero-Coupon', 'Bond', 'CBond', 'EuroCall_UL',
  'EuroPut_UL', 'EuroCall_ZC' or 'EuroPut_ZC'.}

  \item{t}{Date of pricing (has to be an integer)}

  \item{T}{Date of maturity for the option}

  \item{nCoupons}{Number of coupons}

  \item{couponsRate}{Rate of coupons}

  \item{omega}{Recoverables in case of default}

  \item{s}{Date of maturity for the underlying}

  \item{Strike}{Strike for options}
}
\description{
  Get a distribution for any asset price. This method is a
  wrapper for asset specific pricers.
}
\examples{
objScenario  <- new("Scenarios")
# Basic scenario's parameters setting
objScenario  <- setParamsBaseScenarios(objScenario, horizon = 10, nScenarios = 1000)
# Risk factors parameters setting
objScenario  <- setRiskParamsScenariosrt(objScenario, vol = .1, k = 2)
objScenario  <- setRiskParamsScenariosS(objScenario, vol = .1, k = 2,
volStock = .2, stock0 = 100, rho=.5)
objScenario  <- setRiskParamsScenariosliqSpr(objScenario, eta=.05, liquiditySpread0=.01)
objScenario  <- setRiskParamsScenariosdefSpr(objScenario, volDefault=.2,
defaultSpread0=.01, alpha=.1, beta=1)
# Forward and ZC rates setting
data(ZC)
objScenario  <- setForwardRates(objScenario, ZC, horizon=10)
objScenario  <- setZCRates(objScenario, ZC, horizon=10)
# Projection
objScenario  <- customPathsGeneration(objScenario, type="shortRate")
objScenario  <- customPathsGeneration(objScenario, type="stock")
objScenario  <- customPathsGeneration(objScenario, type="defaultSpread")
objScenario  <- customPathsGeneration(objScenario, type="liquiditySpread")
Asset_PriceDistribution(objScenario,type='ConvBond',t=0,T=10,nCoupons=1,couponsRate=0.03)
}

