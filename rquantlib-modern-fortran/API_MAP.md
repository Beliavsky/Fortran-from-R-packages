# API map

## Direct or close numerical mappings

| RQuantLib surface | Modern Fortran procedure | Notes |
|---|---|---|
| `EuropeanOption` | `european_option` | Black-Scholes-Merton value and Greeks |
| `EuropeanOptionArrays` | `european_option_array` | Strike by maturity grid |
| `EuropeanOptionImpliedVolatility` | `european_implied_volatility` | Bracketed inversion |
| `AmericanOption` | `american_option` | CRR tree, not all QuantLib engines |
| `AmericanOptionImpliedVolatility` | `american_implied_volatility` | CRR-based inversion |
| `BinaryOption` | `binary_option` | Cash/asset, European/American |
| `BinaryOptionImpliedVolatility` | Use `european_implied_volatility` with a binary price wrapper | No dedicated wrapper |
| `BarrierOption` | `barrier_option` | Recombining-tree numerical analogue |
| `AsianOption` | `geometric_asian_option`, `arithmetic_asian_mc` | Continuous geometric and MC arithmetic |
| `dayCount` | `day_count` | Four supported bases |
| `yearFraction` | `year_fraction` | Four supported bases |
| `advanceDate`, `advance` | `add_days`, `add_months`, `advance_date` | Plain Gregorian dates |
| Calendar predicates | `weekday`, `is_weekend`, `is_business_day` | Weekend plus custom holidays |
| Calendar adjustment | `adjust_date` | Five conventions |
| `businessDaysBetween` | `business_days_between` | Plain calendar object |
| `Schedule` | `make_schedule` | Regular forward/backward schedules |
| `DiscountCurve` flat input | `make_flat_curve` | Continuous flat zero rate |
| `DiscountCurve` zero quotes | `make_zero_curve` | Log-linear discounts |
| Deposit/future/swap curve construction | `bootstrap_curve` | Simplified aligned bootstrap |
| `ZeroPriceByYield` | `zero_price_by_yield` | Simple/compounded/continuous |
| `ZeroYield` | `zero_yield_by_price` | Closed-form inversion |
| `FixedRateBondPriceByYield` | `fixed_rate_bond_from_yield` | Regular coupon schedule |
| `FixedRateBondYield` | `fixed_rate_bond_yield` | Bracketed inversion |
| `FixedRateBond` | `fixed_rate_bond_from_curve` | Curve-based pricing |
| `FloatingRateBond` | `floating_rate_bond_from_curves` | Forward/discount multicurve analogue |
| `FittedBondCurve` Nelson-Siegel | `fit_nelson_siegel` | Variable projection |
| `FittedBondCurve` Svensson | `fit_svensson` | Variable projection |
| SABR smile calculation | `sabr_lognormal_volatility` | Hagan lognormal formula |
| `SabrSwaption` calibration core | `calibrate_sabr` | Fixed-beta smile calibration, not cube machinery |
| Hull-White cap calibration core | `calibrate_hull_white_caplets` | Plain caplet-price arrays |
| Hull-White affine pricing | `hull_white_discount_bond`, `hull_white_bond_option` | One-factor formulas |
| European swaption core | `hull_white_swaption` | Jamshidian decomposition |

## Fortran-specific result types

- `option_result`
- `discount_curve_t`
- `curve_fit_result`
- `bond_result`
- `sabr_result`
- `hull_white_calibration_result`
- `date_t`, `calendar_t`, and `schedule_t`

These replace R lists, S3 classes, date vectors, and Rcpp return objects.

## Explicitly not mapped

- `BermudanSwaption`
- `AffineSwaption` G2++ and Black-Karasinski variants
- `CallableBond`
- `ConvertibleZeroCouponBond`
- `ConvertibleFixedCouponBond`
- `ConvertibleFloatingCouponBond`
- Full QuantLib calendar catalogue
- Complete piecewise-yield-curve helper and interpolation catalogue
- Full SABR volatility cubes
- Rcpp modules and capability/version queries
- Plot, print, and summary methods
