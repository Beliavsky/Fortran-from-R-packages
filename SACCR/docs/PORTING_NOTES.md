# Porting notes

## Translation approach

The R implementation builds a mutable `data.tree` and then walks that tree to
calculate add-ons. The Fortran implementation calculates the same hierarchy
directly from `trade_t` arrays and stores its output in typed arrays. This
removes a non-computational dependency and makes invalid field names or node
shapes compile-time rather than runtime concerns.

The attached modern Fortran `Trading` translation is bundled as an independent
FPM path dependency. SACCR uses its `trade_t`, `csa_t`, `collateral_t`, CSV
readers, supervisory-duration, maturity-factor, and option-delta methods.

## Documented upstream corrections

The following apparent defects or inconsistent edge cases were corrected
rather than reproduced:

1. **IRD bucket cross-term.** The upstream expression repeats the bucket-1 /
   bucket-2 product in its final `0.6` term. The port uses the standard
   bucket-1 / bucket-3 term:
   `D1^2 + D2^2 + D3^2 + 1.4 D1 D2 + 1.4 D2 D3 + 0.6 D1 D3`.

2. **Commodity netting.** The upstream code applies `abs` after every trade,
   making the result depend on input order. The port nets signed effective
   notionals first and applies the absolute value once.

3. **Commodity OEM double count.** The upstream OEM branch adds both effective
   notional and `Notional * Ei`. The port uses the OEM notional-maturity term
   once.

4. **Simplified commodity correlation.** The upstream routine overwrites the
   supervisory-data correlations for simplified SA-CCR but later hardcodes
   commodity correlation to `0.4`. The port uses correlation `1` in simplified
   mode.

5. **Prefixed FX and other-exposure grouping.** Some upstream grouping filters
   compare original trade fields against already prefixed `Basis_` or `Vol_`
   node names, yielding empty groups. The port keeps the grouping key separate
   from the underlying trade field.

6. **Digital-option cap state.** The upstream conditional initializes and
   resets the exotic identifier in an order that prevents the intended paired
   cap from being applied reliably. The port groups expanded digital legs by
   identifier and applies the cap to their signed aggregate.

7. **Percentage CSA values.** The R documentation permits percentage values,
   but `CalcRC` leaves required variables undefined for that mode. The port
   interprets them as proportions of absolute portfolio MtM.

8. **`ignore_margin`.** The upstream branch can still subtract collateral even
   when margin is requested to be ignored. The port treats the exposure as
   fully unmargined.

9. **Current FX hedge amount.** The upstream example subtracts a
   risk-weighted CDS amount from EAD. The port separates protected and
   unprotected EAD first, then applies the two risk weights.

10. **Zero add-on PFE.** The port returns zero directly when aggregate add-on is
    zero, avoiding a possible `0/0` multiplier.

## Numerical compatibility

The tests reproduce the documented package examples:

| Example | Fortran EAD |
|---|---:|
| IRD | 569.47014094 |
| Credit | 381.23831875 |
| Commodity | 5405.61598246 |
| FX | 924.00000000 |
| IRD plus credit | 936.45050554 |
| Margined IRD plus commodity | 1879.21263150 |

The R documentation rounds these values to 569, 381, 5406, 936, and 1879.

## FPM and compiler scope

The source uses Fortran 2018 syntax, lower-case identifiers, `implicit none`,
and `dp = kind(1.0d0)` inherited from the Trading dependency. There are no
external compiled dependencies beyond the bundled package.
