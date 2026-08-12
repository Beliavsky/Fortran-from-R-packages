# genalg-fortran

Modern Fortran translation of the computational core of the R package
`genalg` 0.2.1 by Egon Willighagen and Michel Ballings.

The original package implements two simple minimization-oriented genetic
algorithms entirely in R:

- `rbga()` for bounded floating-point chromosomes;
- `rbga.bin()` for binary chromosomes.

This project translates those algorithms to Fortran 2018 and packages them
with FPM. R plotting, S3 summary/plot methods, and R-specific object handling
are intentionally omitted.

## Implemented behavior

The translation retains the numerical behavior that defines the package:

- minimization of a user callback;
- random or user-supplied initial populations;
- default mutation probability `1/(nvar+1)`;
- default elitism `floor(pop_size/5)`;
- generation-by-generation best and mean objective histories;
- cached objective values for chromosomes copied unchanged;
- rank-biased parent selection using the same normal-density weights as the
  R code;
- sampling two distinct parents without replacement;
- one-point crossover with crossover points from `0:nvar`, including exact
  parent copies at the endpoints;
- the special one-variable no-crossover branch;
- bounded real mutation with the same generation damping and fallback to a
  fresh uniform value when the mutation leaves the domain;
- binary initialization according to `zeroToOneRatio`;
- rejection of all-zero chromosomes during binary random initialization and
  after mixed-parent crossover;
- binary mutation by resampling a bit with the same zero/one weighting;
- optional monitor callbacks after each evaluated generation.

The real-valued mutation deliberately preserves the original expression

```text
mutationVal = stringMax[var] - stringMin[var] * 0.67
```

including its operator precedence. This is not changed to
`0.67*(stringMax-stringMin)`.

## Binary cached-fitness compatibility

`genalg` 0.2.1 has a small computational quirk in `rbga.bin()`: when mutation
changes a chromosome that inherited an already-known evaluation through an
endpoint crossover, the cached evaluation is not invalidated. The chromosome
can therefore carry a stale fitness into the next generation.

For compatibility, the Fortran default is:

```fortran
control%legacy_binary_eval_cache = .true.
```

For new applications, the corrected behavior is recommended:

```fortran
control%legacy_binary_eval_cache = .false.
```

A regression test exercises both modes explicitly.

## Public API

The main module is `genalg`.

```fortran
use genalg, only : dp, rbga_control, rbga_result, rbga
```

A floating-point objective has the interface

```fortran
function objective(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: f
end function objective
```

and is called with

```fortran
call rbga(lower, upper, objective, result, control)
```

Optional suggestions are supplied as a matrix with one chromosome per row:

```fortran
call rbga(lower, upper, objective, result, control, suggestions)
```

The binary interface is analogous:

```fortran
use genalg, only : rbga_bin_control, rbga_bin_result, rbga_bin
call rbga_bin(nbits, objective, result, control)
```

Binary objective callbacks take `integer :: x(:)` containing zeroes and ones.

Results contain the final population and evaluations, best/mean histories,
a convenience best chromosome/value, the actual number of objective callback
calls, and the total number of gene mutations.

## RNG compatibility

The R package uses R's global RNG. This port uses a small deterministic
Park-Miller RNG so that it is standalone and reproducible without R.
Consequently, the same integer seed does **not** produce the same stochastic
trajectory as R. The algorithmic rules and sampling distributions are kept,
not R's exact random-number stream.

## Build

With FPM:

```text
fpm build
fpm test
fpm run --example pi_sqrt50
fpm run --example binary_selection
```

The project has no external library dependency.

## Validation

The release tests cover:

- real optimization on the package documentation's `(pi, sqrt(50))` target;
- 30-bit binary optimization;
- suggestions and monitor callbacks;
- exact legacy versus corrected binary evaluation-cache behavior;
- one-gene real and binary branches.

The release is also compiled directly with GNU Fortran using

```text
-std=f2018 -O0 -g -fcheck=all -Wall -Wextra \
-Wimplicit-interface -Werror=implicit-interface
```

and again with `-O2`.

## License and provenance

`genalg` 0.2.1 declares `License: GPL-2`. This translation is therefore
released under GPL-2.0-only. The original package metadata, R sources,
manual pages, and changelog are retained under `original/` for provenance.
See `LICENSES.md` and `TRANSLATION_NOTES.md`.
