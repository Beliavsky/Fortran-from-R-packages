# qap-fortran

Modern Fortran/FPM translation of the computational code in the R package
`qap` 0.1-2 by Michael Hahsler and contributors.

The package solves symmetric nonnegative Quadratic Assignment Problems (QAPs)
with the simulated-annealing heuristic of Burkard and Rendl (1984), translated
from the legacy `qapsim.f` routine shipped with the R package.

## Scope

Translated computational functionality:

- `qap.obj` -> `qap_obj`
- `qap(..., method="SA")` -> `qap_solve`
- Burkard-Rendl simulated annealing, including repetitions/restarts
- the symmetric-QAP O(n) swap-delta update
- `read_qaplib`
- QAPLIB problem and solution data bundled with the original package

R registration, `.Fortran` glue, R RNG wrappers, printing, and R attributes are
not needed in the standalone Fortran library and are omitted.

## Build

With FPM:

```sh
fpm build
fpm test
fpm run --example solve_had20
```

The implementation is standard Fortran 2018 and has no external numerical
library dependencies.

## API

```fortran
use qap

type(qap_problem_t) :: p
type(qap_control_t) :: ctl
type(qap_result_t) :: res

call read_qaplib('data/qaplib/had20.dat', p)
ctl%rep = 10
ctl%seed = 1000_i64
call qap_solve(p%A, p%B, res, ctl)

print *, res%objective
print *, res%permutation
```

The default controls mirror the R package:

```text
rep      = 1
miter    = 2*n       (represented by miter=-1 until n is known)
fiter    = 1.1
ft       = 0.5
maxsteps = 50
```

`qap_result_t` additionally reports attempted swaps, accepted swaps, duplicate
index trials, total cooling steps, and which restart produced the best result.

## Numerical fidelity

The translated annealer preserves several details of the original Fortran code
which are easy to lose in a generic rewrite:

- trial indices use `INT(u*n + 0.5)`, including its endpoint bias;
- the iteration count after cooling is truncated back to integer;
- uphill moves with `delta/T > 10` are assigned zero acceptance probability;
- the stopping range is updated only after accepted swaps or equal-index trials;
- the initial incumbent threshold is the original rough mean-objective estimate;
- repetitions start from independent random permutations, as in the R wrapper.

The R package obtains both the initial permutation and annealing uniforms from
R's RNG. The standalone port uses an explicit Park-Miller RNG so `seed` is
portable and reproducible, but it is not intended to reproduce `set.seed()`'s
exact R random sequence.

## QAPLIB reader

`read_qaplib` reads the `.dat` file and automatically reads a sibling `.sln`
file when present. It intentionally returns the solution permutation exactly as
stored, matching the R package. Some historical QAPLIB solution files use the
inverse assignment convention; the reader does not silently invert them.

All 136 `.dat` files bundled with the package were parsed successfully during
release validation.

## Validation

The test suite checks:

- the known `had12` objective 1652;
- the O(n) swap-delta formula against complete objective recomputation;
- deterministic seeding;
- representative QAPLIB formats, including large instances;
- the README `had20` case, where 10 restarts with the release seed reach the
  known optimum 6922.

See `VALIDATION.md` for compiler flags and details.

## License

The original R package declares GPL-3. This translation is distributed under
GPL-3.0-only. The GPL v3 text is in `LICENSE`; the complete attached source tree
is retained under `original/qap-master/` for provenance.

## Reference

R. E. Burkard and F. Rendl (1984), "A thermodynamically motivated simulation
procedure for combinatorial optimization problems," European Journal of
Operational Research 17(2), 169-174.
