# Computational API coverage

Status: substantial. The exported computational denominator is exactly two:
`priceCDS` and `bootstrapCDS`, as listed in upstream/NAMESPACE. Both are mapped
(2 of 2, 100%). No computational export is omitted. Upstream has no additional
exported computational S3 methods. Helpers are not added to the denominator.

| R export | Implementation | Scope |
| --- | --- | --- |
| priceCDS | credule:price_cds | Full yield/survival interpolation, separate premium/default grids, recovery, accrual, multi-maturity pricing; direct linear leg ratio |
| bootstrapCDS | credule:bootstrap_cds | Full sequential interval-hazard bootstrap, original Brent algorithm and [0,30] bracket, survival output |

Both mappings are **substantial**, not full interface compatibility. Fortran
outputs are arrays, not R data frames; tenors remain with the caller. Explicit
integer status codes replace R warnings/stops. Scalar options replace R vector
defaults that upstream accesses only through their first element. No R recycling,
implicit type coercion or R-specific marshaling is provided.

Internal computation coverage:

- getDiscountFactor -> discount_factor / private discount_raw
- getSurvivalProbability -> survival_probability / private survival_raw
- calculatePremiumLeg / calculateDefaultLeg -> private leg_values / public cds_legs
- objfunFindParSpread -> direct protection/premium-annuity ratio
- objfunFindHazardRate -> local evaluations inside bootstrap_cds
- zero -> pure reverse-communication Brent routines in credule_brent

The C allocation helpers and debug printing are replaced by caller-sized arrays
and local state. Only `price_cds` and `bootstrap_cds` count toward export coverage.

Preserved quirks: right-endpoint default discounting, half-period premium accrual,
truncation of both payment grids, zero-hazard pricing beyond the supplied curve,
and separate leg-specific bootstrap anchors at non-grid knots. See README for
the resulting fractional-tenor round-trip caveat.

Deliberate improvements: explicit input and extent checks, quiet-NaN failure
outputs, exact zero/upper-endpoint hazard acceptance, bounded Brent iteration,
nonfinite/underflow detection, rejection of unbracketed roots, and no artificial
par-spread upper bound. Maturity zero, empty curves, survival zero, or a maturity
shorter than either leg interval are invalid rather than relying on undefined C
behavior. Full-recovery pricing is valid, but full-recovery calibration is not
identifiable and is rejected.
