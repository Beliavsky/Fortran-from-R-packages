# DiceDesign-fortran

Modern Fortran/FPM translation of the computational code in **DiceDesign 1.10**
(Designs of Computer Experiments).

The library implements space-filling designs, Latin-hypercube construction and
optimization, distance/discrepancy criteria, randomized-sphere statistics, and
uniformity tests without requiring R at runtime.

## Implemented functionality

- full factorial, Latin hypercube, orthogonal LHS, NOLH and NOLHDR designs;
- Faure sequence and Faure-prime designs;
- scaling, inverse scaling, empirical uniformization, and numeric xDRDN transforms;
- minimum distance, coverage, mesh ratio, phi-p, and minimum-spanning-tree criteria;
- C2, L2, L2-star, M2, S2, W2, and mixed-L2 discrepancies;
- discrepancy and maximin LHS optimization with simulated annealing and ESE;
- D-max, Strauss, and WSP stochastic designs;
- Greenwood, maximum-spacing, Quesenberry-Miller, Kolmogorov-Smirnov, Kuiper/V,
  and Cramer-von Mises uniformity statistics;
- 2-D and 3-D randomized-sphere statistics;
- the complete upstream NOLH/NOLHDR design tables through dimension 29 and the
  Greenwood critical-value table.

Plotting and R-specific formula/list/data-frame/display infrastructure are not
translated. The optional `rgl`, `lattice`, and `randtoolbox` packages are not
runtime dependencies of this port.

## Build with FPM

```text
fpm build
fpm test
fpm run --example basic_designs
```

The translation targets Fortran 2018. No external numerical library is required.

## Minimal example

```fortran
program demo
  use iso_fortran_env, only : int64
  use dicedesign, only : dp, lhs_design, discrepancy_value, &
    discrep_ese_lhs, lhs_optimization_result
  implicit none

  real(dp), allocatable :: x(:, :)
  type(lhs_optimization_result) :: opt

  call lhs_design(12, 3, x, randomized=.false., seed=12345_int64)
  call discrep_ese_lhs(x, opt, inner_iterations=20, candidates=12, &
    outer_iterations=2, criterion='C2', seed=12345_int64)

  print *, discrepancy_value(x, 'C2')
  print *, discrepancy_value(opt%design, 'C2')
end program demo
```

## RNG behavior

The original R functions use R's global RNG. This port uses an explicit local
`rng_state`, so independent reproducible streams are possible without global
state. Seeds therefore do **not** reproduce R's `set.seed()` stream bit-for-bit;
the design algorithms and acceptance rules are preserved.

## License

The upstream package declares `GPL-3`. This translation remains
**GPL-3.0-only**. See `LICENSE`, `UPSTREAM.md`, and the original metadata under
`upstream/`.
