# caRamel-fortran

Modern Fortran translation of the computational core of the R package **caRamel 1.5** (2026-05-07), packaged for the Fortran Package Manager (FPM).

The upstream package implements a multi-objective evolutionary optimizer combining directional search in objective-space simplexes with epsilon-dominance/archive ideas. This port removes the R runtime dependency and exposes the numerical algorithms as ordinary Fortran procedures.

## What is translated

The following exported computational routines have native Fortran counterparts:

- `pareto`, `dominate`, `dominated`
- `val2rank`, `boxes`, `downsize`, `decrease_pop`
- `vol_splx`, `Dimprove`
- `rselect`
- `matvcov`, `Cusecovar`
- `Cinterp`, `Cextrap`, `Crecombination`
- `newXval`
- the main `caRamel` optimization loop, exposed as `caramel_optimize`
- the optional first-order sensitivity calculation

The upstream `geometry::delaunayn()` dependency used by `newXval` is replaced by a native n-dimensional incremental Delaunay triangulator (`caramel_delaunay`). No R, C, C++, Qhull, BLAS, or LAPACK dependency is required.

## Intentionally omitted

The request was to skip plotting code, so these R functions are not ported:

- `plot_caramel`
- `plot_pareto`
- `plot_population`

R-specific infrastructure is also omitted: S3/package startup behavior, progress bars, R cluster/parallel dispatch, R global-variable conventions, and listing-file output. The Fortran callback receives the parameter vector directly, so there is no equivalent of the R package's global matrix named `x`.

## Build

```text
fpm build
fpm test
fpm run --example schaffer
```

The project has no external FPM dependencies.

## Main API

```fortran
use caramel, only: dp, caramel_options, caramel_result, caramel_optimize

type(caramel_options) :: options
type(caramel_result) :: result

options%popsize = 100
options%archsize = 100
options%maxrun = 1000
options%repart_gene = [5, 5, 5, 5]

call caramel_optimize(nobj, nvar, minmax, bounds, prec, objective, result, options)
```

The objective callback has the interface

```fortran
subroutine objective(x, values)
    use caramel, only: dp
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: values(:)
end subroutine objective
```

`minmax(j)=.true.` means maximize objective `j`; `.false.` means minimize it, matching the R package.

Set `options%nout > nobj` when the callback returns additional diagnostic values. They are carried through the archive/result but only the first `nobj` columns participate in optimization.

`options%sensitivity=.true.` enables the upstream-style forward finite-difference Jacobians.

Parameter recombination blocks use the public `index_block` derived type.

## Result type

`caramel_result` contains:

- `success`, `message`
- `parameters`: final Pareto parameter archive
- `objectives`: final objective values (plus optional extra callback outputs)
- `derivatives`: `[front, variable, objective]` finite-difference derivatives when requested
- `save_crit`: per-generation call count and best value of each objective
- `total_pop`: final archive + retained non-front population
- `gpp`, `nrun`, `ngen`

As in the R code, `maxrun` is a generation boundary rather than a hard truncation: a final generation may make the call count exceed `maxrun`.

## Numerical/porting notes

See `ALGORITHM_NOTES.md` for details of the Delaunay implementation and documented corrections made while preserving the intended caRamel algorithm. `API_MAPPING.md` maps the R API to the Fortran procedures, and `VALIDATION.md` records strict-build and randomized differential tests.

## License and provenance

The upstream `DESCRIPTION` declares `License: GPL-3 | file LICENSE`; the supplied upstream `LICENSE` file contains the GNU Lesser General Public License version 3. Both are preserved here verbatim (`DESCRIPTION.upstream` and `LICENSE`), and a copy of GPLv3 is included as `LICENSE-GPL-3.txt` for completeness.

This is a translation/derived implementation of caRamel. See `UPSTREAM.md` for authorship and provenance.
