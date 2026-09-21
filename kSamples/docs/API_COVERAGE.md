# API coverage

Coverage is measured over distinct exported computational R functions in `upstream/NAMESPACE`. Plotting, printing, formatting, S3 presentation, and interface-only helpers are excluded. The machine-readable source of truth is `[extra.translation]` in `fpm.toml`.

All 14 exported computational functions have a meaningful Fortran mapping. This is **100% mapping coverage**, not a claim of 100% R compatibility. In particular, formula/list dispatch, NA removal and warnings, R object construction, R RNG stream identity, and presentation methods are outside the Fortran API.

Two numerical/interface differences remain material:

1. `ad.pval`: the upstream R routine applies `smooth.spline` interpolation to the tabulated Scholz-Stephens calibration. The Fortran routine preserves the complete tables and uses deterministic linear interpolation in `1/sqrt(m)` followed by linear interpolation in logit-probability space.
2. `SteelConfInt`: asymptotic simultaneous confidence-bound inversion is translated, but the upstream exact/simulated confidence-level refinement is not yet implemented. Non-asymptotic method requests currently return the asymptotic calibration and are marked partial in `fpm.toml`.

Exact/randomization test enumeration is implemented for `ad.test`, `qn.test`, `Steel.test`, `jt.test`, and contingency tables, with exact block combinations where the R package defines them. Simulation paths use the Fortran intrinsic pseudorandom generator and therefore do not reproduce R RNG streams bit-for-bit.
