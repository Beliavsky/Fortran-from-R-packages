# Translation notes

## Scope

`elliptic` 1.5-1 is a pure-R package (no native C/C++ backend).  The translation ports
the mathematical/computational routines to modern Fortran and omits presentation or
R-runtime infrastructure.

### Direct/equivalent computational mappings

| R routine/family | Fortran |
|---|---|
| `K.fun`, `nome`, `nome.k` | `k_complete`, `nome`, `nome_k` |
| `theta1`--`theta4`, `theta1dash*` | `theta1_q/m`--`theta4_q/m`, `theta1dash*_q/m` |
| `theta.s/c/d/n` | `theta_s/c/d/n` |
| `sn`, `cn`, `dn`, quotient family | same short Fortran names |
| `e16.28.*` | `e16_28_*` |
| `e16.36.*` | implemented by `theta_s/c/d/n` |
| `e16.37.*`, `e16.38.*` | `e16_37_*`, `e16_38_*` |
| `parameters` | `elliptic_parameters`, `parameters_from_g/omega` |
| `e1e2e3`, `eee.cardano` | `cubic_roots`, `eee_cardano` |
| `half.periods`, `as.primitive`, `is.primitive` | `half_periods`, `primitive_periods`, `is_primitive` |
| `mn`, `fpp` | `mn_coordinates`, `fundamental_parallelogram` |
| `P`, `Pdash` | `wp`, `wp_prime` |
| `sigma`, `zeta` | `weierstrass_sigma`, `weierstrass_zeta` |
| `coqueraux` | `coqueraux` |
| `ck`, `amn` | `ck_coefficients`, `amn_coefficients` |
| Laurent routines | `wp_laurent`, `wp_prime_laurent`, `sigma_laurent`, `sigma_prime_laurent`, `zeta_laurent` |
| `g2.fun`, `g3.fun`, `g.fun` | `g2_from_periods`, `g3_from_periods`, `g_from_periods` |
| direct/divisor/fixed/Lambert g2/g3 variants | `g2_direct/divisor/fixed/lambert`, corresponding g3 routines |
| `eta`, `eta.series` | `dedekind_eta`, `eta_series` |
| `J`, `lambda` | `j_invariant`, `modular_lambda` |
| `mob`, `%mob%` | `mobius_transform` |
| `farey`, `unimodular` | `farey_series`, `unimodular_matrices` |
| arithmetic helpers | `primes_upto`, `factorize`, `divisor_sigma`, `totient`, `mobius`, `liouville` |
| `myintegrate`, contour/segment integration, `residue` | `integrate_real_path`, `integrate_segments`, `residue` |
| `newton_raphson` | `newton_raphson_complex` |
| special parameter sets | `equianharmonic_parameters`, `lemniscatic_parameters`, `pseudolemniscatic_parameters` |

## Numerical implementation choices

- The complete elliptic integral uses the same AGM construction as the R code, but the
  Fortran convergence test is based on a relative complex magnitude instead of
  `all.equal()`.
- Theta sums use a term-size convergence test with a large default iteration ceiling.
  This is important for nomes close to the unit circle and covers the upstream issue-7
  regression cases.
- The cubic roots of `4*x^3-g2*x-g3` are obtained by a complex Durand-Kerner solve and
  then ordered using the same period/Coqueraux logic used by the R package.  The
  explicit upstream Cardano formula is separately available as `eee_cardano`.
- `parameters_from_omega()` preserves the supplied primitive period basis; it does not
  round-trip through invariants and silently choose another equivalent basis.
- Contour integration uses complex Simpson integration directly.  R's implementation
  integrated real and imaginary parts separately via `stats::integrate()`.
- `farey_series()` constructs reduced fractions directly, removing the upstream MASS
  dependency used only for rational/fraction formatting.

## Deliberately omitted R/external-system code

- `view()` and `latplot()` are plotting/presentation routines.
- `P.pari()` shells out to the separately installed `pari/gp` executable.  It is an
  external verification bridge rather than a numerical implementation in `elliptic`.
- `Re<-`, `Im<-`, class tags, attributes, `massage()`, and vector-shape restoration are
  R object-system conveniences.  Fortran uses native `complex(dp)` values and an
  explicit `elliptic_parameters` type.
- `limit()` is a display-oriented clipping helper used with visualizations.
- Vectorized R wrappers map naturally to caller loops/array expressions in Fortran;
  the scalar numerical kernels are the canonical API.

## Validation

The tests include fixed upstream reference values for Jacobi functions, Weierstrass P
and sigma, the page-656 Laurent coefficients, page-668 sigma expansions, modular
identities, invariant/period round trips, and contour/residue calculations.
