# Validation

Validated with GNU Fortran using:

```text
-std=f2008
-Wall -Wextra
-Wimplicit-interface -Werror=implicit-interface
-fcheck=all
-ffpe-trap=invalid,zero,overflow
-O0 -g
```

The mstate test suite covers:

- transition topology, `trans2tra`/`tra2trans`, and upstream-compatible path prefixes
- the canonical six-subject `msprep` example
- `crprep` competing-event expansion, stratum-specific censoring weights, truncation, and keep columns
- event counts, cross-sections, landmark cutting, and covariate expansion
- `msdata`/`etm` round trips and cluster bootstrap
- Breslow and Efron tied-event cumulative hazards and variances
- Cox-array and supplied-survival Cox bridges
- reduced-rank alternating Cox fitting and rank-factorization invariants
- Markov-test no-covariate and Cox-adjusted paths
- forward Aalen-Johansen probabilities and restricted ELOS
- single- and multi-start landmark Aalen-Johansen variance
- cumulative incidence point estimates, grouped fits, and influence-function SEs
- forward/reset multi-state sampling, history-state effects, censoring, path and data outputs
- relative-survival transition/hazard/covariance splitting
- Markov-test weight helper dimensions and non-negativity

The complete user-supplied survival source bundle was also compiled separately
under the strict runtime flags in the preceding integration pass. Its supplied
AFT, AJ, counting-process Cox, Cox, Kaplan-Meier, penalized-spline, statistics,
utility, and spline tests passed.
