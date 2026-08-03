# API map

| R export | Fortran counterpart | Coverage |
|---|---|---|
| `coupon.dates` | `coupon_dates`, `coupon_dates_from_years` | Direct typed counterpart, including regular and short/long first/last schedules |
| `coupons` | `coupons`, `coupon_cashflows` | Direct cash-flow calculation with day-count fractions and redemption |
| `discount.factors` | `discount_factors`, `discount_factor` | Continuous and periodically compounded rates |
| `discount.time` | `discount_time` | Direct Quantil-style full-year plus fractional-year calculation |
| `accrued.interests` | `accrued_interests`, `accrued_interest` | High- and low-level counterparts |
| `spot2forward` | `spot2forward`, `spot_to_forward` | Constant and linear conventions |
| `fwd2spot` | `fwd2spot`, `forward_to_spot` | Constant and linear integration |
| `valuation.bonds` | `valuation_bonds` | Scalar/vector rate and full curve overloads |
| `valuation.swaps` | `valuation_swaps`, `swap_leg_value` | Fixed/floating IRS and cross-currency leg combinations |
| `price.dirty2clean` | `price_dirty2clean`, `dirty_to_clean` | Direct counterpart |
| `bond.price2rate` | `bond_price2rate`, `bond_price_to_rate` | Bounded bisection yield inversion |
| `sens.bonds` | `sens_bonds`, `bond_sensitivity` | Modified duration plus richer convexity/DV01 result |
| `average.life` | `average_life` | Cash-flow or discounted weighted average life |
| `curve.calibration` | `curve_calibration`, `bootstrap_curve` | Yield interpolation and price-based constrained calibration separated explicitly |
| `curve.calculation` | `curve_calculation` | Multi-date/multi-column curve construction |
| `basis.curve` | `basis_curve` | Constrained basis-node calibration from swap market values |

## Additional utilities

- ISO-style date parsing and formatting
- Month-end-safe date arithmetic
- weekend business-day conventions `F`, `B`, `MF`, and `MB`
- `30/360`, `ACT/360`, `ACT/365`, `ACT/365L`, `NL/365`, `ACT/ACT-ISDA`, and `ACT/ACT-AFB`
- constant/linear curve interpolation with flat extrapolation
- native bounded Nelder-Mead optimization
