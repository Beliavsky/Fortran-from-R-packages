# Translation notes

## Design

The translation exposes numerical kernels directly through modern Fortran
modules and explicit arrays. The umbrella module is `actuar`. R S3/formula
infrastructure is intentionally omitted.

## v0.3 algorithm provenance

- Minimum-distance objectives follow `R/mde.R`. The Fortran release uses a
  bounded Nelder-Mead optimizer for all parameter dimensions rather than
  reproducing R's `optim` method-selection machinery.
- Coverage transformations follow `R/coverage.R`, including endpoint masses,
  franchise deductibles, per-loss/per-payment conditioning and Jacobian
  scaling under inflation/coinsurance.
- Hachemeister barycenter fitting follows `R/hache.barycenter.R`. Weighted QR
  is implemented with modified Gram-Schmidt; the returned transition matrix
  maps the orthogonalized regression basis back to the original design basis.
- `hierarc_exact_fit` follows both `R/hierarc.R` and `src/hierarc.c`. The
  compiled iteration's fallback from zero credibility to current node weight
  is preserved exactly.
- `rcomphierarc_simulate` follows the hierarchy-expansion and ancestor-parameter
  propagation of `R/rcomphierarc.R`. R expressions are replaced by typed
  procedure callbacks.

The relevant original files are retained under `upstream/reference-v03`.

## Earlier algorithm provenance

The following direct ports remain from v0.1-v0.2:

- Poisson-inverse-Gaussian recurrence from `src/poisinvgauss.c`;
- Panjer recursion from `src/panjer.c` and `R/panjer.R`;
- phase-type formulas from `src/phtype.c`;
- incomplete-beta continuation from `src/betaint.c`;
- exact aggregate convolution from `R/exact.R`;
- normal-power aggregation from `R/normal.R`;
- phase-type ruin fixed point from `R/ruin.R`;
- grouped-data kernels from `R/emm.R`, `R/ogive.R` and `R/elev.R`.

## expint dependency

Upstream `actuar` imports `expint`. The translated `expint-fortran` package is
vendored as an FPM path dependency and supplies extended incomplete-gamma
semantics needed by inverse-family limited moments.

## Threading note for MDE

The MDE optimizer stores its active callback/objective context at module scope
to avoid compiler trampolines from nested procedure callbacks. Consequently,
concurrent calls to the MDE routines from multiple threads require external
serialization in v0.3.
