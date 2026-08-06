# Porting notes

## Scope

All exported computational routines in xVA 1.3 are represented. The package's
IRS simulation path remains restricted to `IRDSwap`, matching the original
package documentation. All CSV data files are retained. `load_supervisory_cva_data` materializes all
eleven supervisory tables, including the sector matrices and commodity, equity,
and hedge-correlation records that are not consumed by the current capital
formula implementation.

## Deliberate corrections

Several upstream expressions appear to be implementation defects rather than
model choices:

1. `CalcNGR` used `sum(max(MtM_Vector, 0))`, where R's scalar `max` returns only
   the largest positive trade. The Fortran routine uses the standard gross
   positive exposure, `sum(max(MtM_Vector, 0))` element by element.
2. The CEM add-on summed percentage factors without multiplying by trade
   notionals. The Fortran calculation uses `abs(notional) * add-on factor`.
3. Duplicate time-grid points are removed. Keeping duplicates can create zero
   time steps and divisions by zero in the forward-rate calculation.
4. SA-CVA interest-rate bumps use common random numbers. The base and bumped
   exposure simulations therefore differ because of the curve bump, not an
   unrelated Monte Carlo draw. The end-to-end calculator uses seed `104729`
   unless the caller supplies another seed.
5. Under IMM, each trade's updated `mtm` is its own first simulated path value.
   The R loop assigned the cumulative portfolio value to every successive
   trade.
6. CEM notionals and non-IMM effective-maturity weights use absolute notionals,
   avoiding cancellation between signed positions.

These changes are exercised by the tests and keep the public mathematical
intent of the package while avoiding unstable or dimensionally inconsistent
results.

## Preserved behavior

The exposure summaries preserve the R package's conditional definitions:

- EE is the mean over nonnegative simulated values;
- NEE is the mean over negative simulated values and remains negative;
- PFE is the type-7 quantile over strictly positive values.

These differ from unconditional expected positive/negative exposure
conventions but are retained for compatibility.

## Regulatory dependencies

SA-CCR calculations call the bundled SACCR Fortran package. Curves, trades,
CSAs, collateral, interpolation, random-number utilities, and type-7 quantiles
come from the bundled Trading Fortran package. The dependency packages remain
separate FPM packages and retain their original licenses.
