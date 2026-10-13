# API coverage and exclusions

Status: **substantial**. **20 of 20 (100%)** exported computational R functions
are mapped. The fraction measures mapped computational functions, not complete
R compatibility, corrected finance theory, or identical stochastic streams.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| setForwardRates | esg_models | esg_forward_rates | substantial |
| rShortRate | esg_models | esg_short_rate | substantial |
| rStock | esg_models | esg_stock | substantial |
| rRealEstate | esg_models | esg_real_estate | substantial |
| rLiquiditySpread | esg_models | esg_liquidity | substantial |
| rDefaultSpread | esg_models | esg_default | substantial |
| rAllRisksFactors | esg_models | esg_generate | substantial |
| rAssetDistribution | esg_pricing | esg_simulate_asset | substantial |
| MartingaleTest | esg_models | esg_martingale | substantial |
| customPathsGeneration | esg_models | esg_generate | substantial |
| Asset_PriceDistribution | esg_pricing | esg_asset_price | substantial |
| Bond_PriceDistribution | esg_pricing | esg_bond | substantial |
| CDSPremium_PriceDistribution | esg_pricing | esg_cds | substantial |
| ConvBond_PriceDistribution | esg_pricing | esg_convertible_bond | substantial |
| CBond_PriceDistribution | esg_pricing | esg_corporate_bond | substantial |
| EuroCall_Stock_PriceDistribution | esg_pricing | esg_stock_call | substantial |
| EuroPut_Stock_PriceDistribution | esg_pricing | esg_stock_put | substantial |
| EuroCall_ZC_PriceDistribution | esg_pricing | esg_zc_call | substantial |
| EuroPut_ZC_PriceDistribution | esg_pricing | esg_zc_put | substantial |
| ZCBond_PriceDistribution | esg_pricing | esg_zc_bond | substantial |

## Mapping limitations

- **setForwardRates**: Finite-difference extraction retained; typed arrays replace S4 setter and strict shape checks replace R NA removal.
- **rShortRate**: Full short-rate recurrence with explicit innovations and extracted forwards; R object construction and RNG omitted.
- **rStock**: Full correlated stock recurrence using rate paths and explicit normals; redundant discarded R rate draws are not consumed.
- **rRealEstate**: Full estate recurrence using rate paths and explicit normals; redundant discarded R rate draws are not consumed.
- **rLiquiditySpread**: Full arithmetic liquidity process with explicit innovations rather than R RNG.
- **rDefaultSpread**: Upstream volDefault/4 drift retained; negative paths return status 3 rather than subsequent sqrt NaN warnings.
- **rAllRisksFactors**: All five risks generated into typed scenario arrays with controlled innovations; no R RNG sequence parity.
- **rAssetDistribution**: Required risks generated and all nine asset branches valued; supplied innovations and explicit horizon replace hidden RNG/S4; broken stock-call alias repaired.
- **MartingaleTest**: Mean discounted stocks and maximum deviation computed; plot, lines and title omitted.
- **customPathsGeneration**: All five processes and risk selection supported; typed scenario replaces S4, supplied normals replace rnorm.
- **Asset_PriceDistribution**: All nine pricing branches; typed arguments and status replace R missing/recycling/warnings. Both stock-call aliases accepted.
- **Bond_PriceDistribution**: Full coupon schedule and source accrued-interest convention; integer coupon frequency and consistent terminal shape.
- **CDSPremium_PriceDistribution**: Source 100-step numerator/denominator and omega multiplier retained; stable LMN algebra extends finite range.
- **ConvBond_PriceDistribution**: Source bond plus unit-strike stock call; source terminal-one convention retained in scenario-sized array.
- **CBond_PriceDistribution**: Full source coupon, bare A(T) principal and 100-step recovery formula; stable equivalent LMN algebra extends finite range.
- **EuroCall_Stock_PriceDistribution**: Source lagged stock index, absolute Tau(T) and terminal-one convention retained; zero variance extended by continuous limit.
- **EuroPut_Stock_PriceDistribution**: Source lagged stock index, absolute Tau(T) and terminal-one convention retained; zero variance extended by continuous limit.
- **EuroCall_ZC_PriceDistribution**: Source d1 logarithm placement and absolute H retained; zero H extended by continuous intrinsic limit.
- **EuroPut_ZC_PriceDistribution**: Source d1 logarithm placement and absolute H retained; zero H extended by continuous intrinsic limit.
- **ZCBond_PriceDistribution**: Upstream monthly indexing, exact grid comparisons and discount interpolation preserved; explicit validation and consistent result shapes.

## Exact excluded exported interfaces

The following S4 slot getters/setters merely store or retrieve values and do not
enter the computational denominator. `setZCRates` slices and stores the supplied
curve; it is available as the uncounted `esg_zero_rates` helper. In contrast,
`setForwardRates` performs finite differences and is included. `MartingaleTest`
has a computational return value and is included despite its omitted plot.

- `setParamsBaseScenarios`: S4 parameter/slot accessor or array-slice setter; no numerical algorithm
- `setRiskParamsScenarios`: S4 parameter/slot accessor or array-slice setter; no numerical algorithm
- `setRiskParamsScenariosrt`: S4 parameter/slot accessor or array-slice setter; no numerical algorithm
- `setRiskParamsScenariosS`: S4 parameter/slot accessor or array-slice setter; no numerical algorithm
- `setRiskParamsScenariosRE`: S4 parameter/slot accessor or array-slice setter; no numerical algorithm
- `setRiskParamsScenariosliqSpr`: S4 parameter/slot accessor or array-slice setter; no numerical algorithm
- `setRiskParamsScenariosdefSpr`: S4 parameter/slot accessor or array-slice setter; no numerical algorithm
- `setZCRates`: S4 parameter/slot accessor or array-slice setter; no numerical algorithm
- `getParamsBaseScenarios`: S4 parameter/slot accessor or array-slice setter; no numerical algorithm
- `getRiskParamsScenarios`: S4 parameter/slot accessor or array-slice setter; no numerical algorithm
- `getRiskParamsScenariosrt`: S4 parameter/slot accessor or array-slice setter; no numerical algorithm
- `getRiskParamsScenariosS`: S4 parameter/slot accessor or array-slice setter; no numerical algorithm
- `getRiskParamsScenariosRE`: S4 parameter/slot accessor or array-slice setter; no numerical algorithm
- `getRiskParamsScenariosliqSpr`: S4 parameter/slot accessor or array-slice setter; no numerical algorithm
- `getRiskParamsScenariosdefSpr`: S4 parameter/slot accessor or array-slice setter; no numerical algorithm
- `getForwardRates`: S4 parameter/slot accessor or array-slice setter; no numerical algorithm
- `getZCRates`: S4 parameter/slot accessor or array-slice setter; no numerical algorithm
- `getShortRatePaths`: S4 parameter/slot accessor or array-slice setter; no numerical algorithm
- `getstockPaths`: S4 parameter/slot accessor or array-slice setter; no numerical algorithm
- `getrealEstatePaths`: S4 parameter/slot accessor or array-slice setter; no numerical algorithm
- `getLiquiditySpreadPaths`: S4 parameter/slot accessor or array-slice setter; no numerical algorithm
- `getdefaultSpreadPaths`: S4 parameter/slot accessor or array-slice setter; no numerical algorithm
- `L`: unexported internal helper
- `K`: unexported internal helper
- `Tau`: unexported internal helper
- `A`: unexported internal helper
- `B`: unexported internal helper
- `C`: unexported internal helper
- `G`: unexported internal helper
- `H`: unexported internal helper
- `Beta`: unexported internal helper
- `Alpha`: unexported internal helper
- `CouponTimeValue`: unexported internal helper
- `ForwardExtraction`: unexported internal helper
- `ZCExtraction`: unexported internal helper

S4 class exports ParamsScenarios and Scenarios are represented by typed Fortran
records and are not functions. Internal helpers K/L/Tau/A/B/C/G/H/Alpha/Beta,
CouponTimeValue, ForwardExtraction and ZCExtraction are not separately counted.
No computational export is omitted from the denominator or marked complete.
