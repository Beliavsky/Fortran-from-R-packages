# ldhmm-fortran

A modern Fortran/FPM translation of the computational code in the R package
`ldhmm` 0.6.1, Hidden Markov Models for financial time series with symmetric
lambda (exponential-power/generalized-normal) state distributions.

The project is self-contained and has no external numerical-library dependency.
It retains the upstream package under `original/ldhmm-master/` and preserves the
upstream Artistic-2.0 license.

## Implemented

- Symmetric lambda distribution: density, CDF, CCDF, moments, and random draws
- HMM construction, validation, transition-matrix initialization, and stationary probabilities
- Natural/working parameter transformations
- Minus log-likelihood and log-space forward/backward recursions
- Local posterior decoding and global Viterbi decoding
- Conditional distributions and pseudo-residuals
- State, density, and volatility forecasts
- Theoretical and decoded state statistics
- State/observation simulation and simulated absolute-return ACF
- Simple moving averages, outlier removal, absolute-return ACF, and price-to-log-return conversion
- Maximum-likelihood fitting with finite-difference BFGS and Nelder-Mead

## Deliberately omitted

Plotting, `ggplot2`/graphics integration, FRED downloading, `xts`/`zoo` class
machinery, YAML configuration, RData loading, and package-specific archived data
access are not translated. The archived upstream files remain in `original/` for
reference and license preservation.

## Build

```sh
fpm build
fpm test
fpm run
```

A direct `gfortran` test script is also supplied:

```sh
./scripts/test_gfortran.sh
```

On Windows with `gfortran` in `PATH`:

```bat
scripts\test_gfortran.bat
```

## Minimal example

```fortran
program example
   use ldhmm
   implicit none

   type(ldhmm_model) :: model, decoded
   real(dp) :: param(2,3), gamma_matrix(2,2), delta(2), x(6)
   integer :: status

   param(1,:) = [0.002_dp, 0.015_dp, 1.0_dp]
   param(2,:) = [-0.003_dp, 0.040_dp, 1.4_dp]
   gamma_matrix(1,:) = [0.97_dp, 0.03_dp]
   gamma_matrix(2,:) = [0.08_dp, 0.92_dp]
   delta = [0.70_dp, 0.30_dp]
   x = [0.01_dp, -0.02_dp, 0.005_dp, 0.03_dp, -0.01_dp, 0.02_dp]

   model = ldhmm_create(2, param, gamma_matrix, delta, status=status)
   decoded = ldhmm_decode(model, x, status=status)
   print '(a,*(f9.5,1x))', 'last state probabilities: ', decoded%states_prob(:,size(x))
end program example
```

See `API.md`, `PORTING.md`, and the programs under `example/` for more detail.

## License

Artistic-2.0. See `LICENSE` and `NOTICE.md`.
