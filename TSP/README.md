# TSP-fortran

Modern Fortran/FPM translation of the computational portion of the R package
**TSP 1.2.7** by Michael Hahsler and Kurt Hornik.

## Scope

Implemented natively in Fortran:

- tour representation (`tsp_tour`) and tour-length computation
- symmetric TSP and asymmetric TSP distance matrices
- Euclidean coordinate distance construction (ETSP)
- identity and random tours
- nearest-neighbor and repetitive nearest-neighbor heuristics
- nearest, farthest, cheapest, and arbitrary insertion heuristics
- insertion-cost kernel from the upstream C implementation
- asymmetric-capable 2-opt local search and the upstream symmetric-specialized 2-opt kernel
- simulated annealing with reversal, swap, or mixed local moves
- `solve_tsp` dispatcher, repetitions, optional 2-opt refinement, and controls
- replacement of positive/negative infinities using the upstream penalty rule
- dummy-city insertion for path/clustering transformations
- Jonker-Volgenant ATSP-to-symmetric-TSP reformulation and dummy filtering
- single- and multiple-cut tour-to-path operations
- TSPLIB read/write support for the formats handled by upstream
- TSPLIB `ATT` and `GEO` distance calculations

Intentionally not translated:

- plotting/image methods
- R S3 object dispatch, printing, formula/data-frame conveniences, and labels
- parallel `foreach` plumbing
- the Concorde and Chained-Lin-Kernighan executable wrappers. The solver code
  is not part of the upstream R package; only command-line interfaces are.

The native simulated-annealing implementation follows the upstream controls and
the logarithmic SANN-style temperature schedule, but it is not intended to
reproduce R's RNG stream or `stats::optim()` bit-for-bit.

## 2-opt upstream bug fix

The supplied upstream `src/two_opt.c` omits the old last-to-first edge from the
improvement calculation when a candidate reversal reaches the final city. This
can cause a valid improving move to be missed. The Fortran translation corrects
that closing-edge term. Randomized asymmetric local-optimality tests are
included to protect the fix.

## Build with FPM

```text
fpm build
fpm test
fpm run --example basic
```

The public umbrella module is `tsp`:

```fortran
program example
    use tsp
    implicit none

    real(dp) :: cost(5,5)
    type(tsp_tour) :: tour

    ! Fill a 5 x 5 distance matrix here.
    cost = 0.0_dp

    tour = solve_tsp(cost)
    print *, tour%order
    print *, tour%length
end program example
```

The default solver matches the R package's default strategy:
`arbitrary_insertion` followed by `two_opt`.

## Direct compiler validation

The translation was validated with GNU Fortran 14.2 using, among other flags,
`-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all`.
See `VALIDATION.md` for details.

## License

GPL-3.0-only. See `LICENSE`, `COPYING`, and `UPSTREAM.md`.
