# Validation

The release source was compiled with gfortran 14.2 using:

```
-std=f2018 -O2 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

All eight permanent test programs pass:

1. `test_core` - parameter repair, EI formula, nondomination, exact
   hypervolume, epsilon indicator and ParEGO scalarization.
2. `test_criteria` - mean, SE, EI, CB, AEI, EQI and adaptive-CB all produce
   finite predictions from a fitted GP.
3. `test_single_mbo` - end-to-end EI optimization of a nonlinear 1D target.
4. `test_batch` - three-point constant-liar batches, including duplicate
   filtering and temporary GP updates.
5. `test_mspot_cb` - Exp(1) parallel-CB proposals and native MSPOT-style
   multiobjective proposal.
6. `test_multiobjective` - DIB/SMS optimization and nondominated final front.
7. `test_parego` - multi-point ParEGO with separate scalar surrogate fits.
8. `test_mixed_continue` - real/integer/categorical parameters and continuation.

The 1D optimization test locates its target near `x=0.379` with objective
approximately `9.56e-3` under the fixed validation seed.

## Independent hypervolume differential test

250 random cases were generated independently in Python with 1-6 points and
2-4 objectives. The Fortran recursive hypervolume was compared with an
independent inclusion-exclusion calculation over all subsets of rectangles.
There were zero failures at an absolute tolerance of `2e-12`; the maximum
absolute discrepancy was approximately `5.33e-15`.

## Build notes

The `fpm` executable was not installed in the validation container. The exact
FPM source/dependency/test tree was therefore compiled directly with gfortran.
`fpm.toml` is parsed separately during the release audit. The manifest keeps
implicit typing and implicit external procedures disabled.

## Independent indicator differential test

120 randomized multi-objective cases (2-4 objectives, 1-5 front points and
1-4 candidates) were compared with an independent Python transcription of
the upstream `c_eps_indicator` and `c_sms_indicator` formulas. The additive
epsilon indicator matched exactly at parsed double precision; the maximum
SMS discrepancy was approximately `1.11e-15`. There were zero failures at an
absolute tolerance of `2e-12`.
