# Porting notes

## Object model

R S4 classes were replaced by typed Fortran derived types:

- `arma_copula`, `sarma_copula`
- `pair_copula`
- `dvine_copula`, `dvine2_copula`, `dvine3_copula`
- `tscopula_spec`, `vtscopula_spec`, `tscm_spec`
- `margin_spec` and typed fit-result structures

R formulas, slots, named lists, `xts`/`zoo` indexes, plotting, and formatted S4
methods are not reproduced.

## Numerical replacements

- `FKF::fkf` is replaced by an exact stationary scalar-observation Kalman
  filter for the upstream ARMA state-space representation. The stationary
  covariance is obtained by Lyapunov fixed-point iteration.
- `ltsa::tacvfARMA` is replaced by a converged impulse-response sum.
- `polyroot` is not needed in the compiled library. AR stability is checked by
  the inverse Levinson/reflection-coefficient recursion.
- `rvinecopulib` is replaced by self-contained pair-copula CDF, density,
  h-function, inverse-h, rotation, D-vine, and Rosenblatt algorithms.
- Gaussian and Student pair-copula CDFs use adaptive one-dimensional
  integration. Archimedean pair densities use stable h-function derivatives
  where a compact analytic density was not retained.
- `arfima` and fractional Brownian Kendall-PACF calculations use deterministic
  truncation of their infinite representations.
- `kdensity` is replaced by a Gaussian KDE for the empirical-margin prediction
  surface.

## Fitting

- ARMA and parametric margins use bounded Nelder-Mead optimization.
- Finite D-vines use sequential conditional Kendall estimation rather than the
  upstream joint `optim` fit. This preserves the pair-family structure and is
  deterministic, but estimates need not be identical to R.
- `fit_tscm_steps` implements the upstream IFM/stepwise workflow.
- `fit_full` currently aliases the stepwise workflow; it does not perform the
  upstream single joint optimization over every margin, V-transform, W-copula,
  and core-copula parameter.
- W-copula branch likelihood contributions follow the four-zone formula in the
  upstream `wobjective` implementation. W-copula simulation uses a latent
  first-order pair-copula chain.

## D-vine families

Implemented families are independence, Gaussian, Student-t, Clayton, Gumbel,
Frank, Joe, and BB1, including rotations. Other `rvinecopulib` families are not
included.

## Licensing

The translated library is GPL-3.0-only, matching upstream. Earlier Fortran
translations of imported packages are retained under `provenance/` only and
are not compiled or linked. In particular, the GPL-2-only `polynom` translation
is not combined with this GPL-3-only build; the required numerical operations
were implemented independently.
