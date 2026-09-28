# mbend

Modern free-form Fortran translation of the computational code in the R package
**mbend 1.3.1**, by Mohammad Ali Nilforooshan. The package bends symmetric
matrices with negative eigenvalues to positive-definite form using the weighted
or unweighted methods of Jorjani et al. (2003) and Schaeffer (2014).

The public Fortran API is in `mbend_api`. The main procedure is:

```fortran
call bend(inmat, result, info, wtmat, reciprocal, max_iter, small_positive, method)
```

`result` is a `bend_result` containing the bent matrix, initial/final
eigenvalues, iteration/convergence metadata, and the ordinary and weighted
performance diagnostics returned by the R implementation. Optional arguments
have the same numerical defaults as the R function: `reciprocal=.false.`,
`max_iter=10000`, `small_positive=1e-4`, and `method="hj"`.

## Dependencies

The translation reuses the repository's shared `rfortran-linalg` package for
symmetric eigendecompositions:

```toml
[dependencies]
rfortran-linalg = { path = "../rfortran-linalg" }
```

No BLAS, LAPACK, or ARPACK source is vendored here, and the package does not
request system `-lblas` or `-llapack` links.

When this directory and `rfortran-linalg` are sibling top-level directories,
run:

```text
fpm build
fpm test
fpm run --example mbend_example
```

## Compatibility notes

The numerical `hj` and `lrs` algorithms, optional weights, reciprocal weights,
zero weights for fixed elements, correlation-matrix diagonal handling, and all
reported numerical diagnostics are translated. R matrix row/column names,
messages/warnings, and list construction are not reproduced. Numerical failures
and invalid arguments use integer status codes; nonconvergence is represented by
`result%converged = .false.` while still returning the last iterate, matching the
R routine's warning-and-return behavior.

As in upstream 1.3.1, a matrix is considered a correlation matrix only when all
diagonal elements are exactly one. The upstream R code only bends matrices whose
smallest computed eigenvalue is negative; already-positive-semidefinite or
positive-definite input is reported with `mbend_err_already_pd`.

An all-zero effective weight matrix is rejected with `mbend_err_zero_weights`.
The upstream R code would divide that matrix by zero and does not provide a
useful numerical result.

## Translation coverage

Package status: **substantial**.

**1 of 1 (100%)** exported computational R functions are mapped. This fraction
measures mapped computational functions, not complete compatibility with R
interfaces or object semantics.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| `bend` | `mbend_api` | `bend` | substantial |

See `fpm.toml` for the machine-readable mapping and `API_COVERAGE.md` for scope
and compatibility details.

## Provenance and license

The upstream package metadata declares `License: GPL-3`, while its repository
`LICENSE` file contains the GNU Lesser General Public License version 3. Both
facts are preserved in this translation; see `NOTICE.md` and `upstream/`.
Translated Fortran source is marked `GPL-3.0-only`, following the license field
in the upstream package metadata. The original upstream license file is retained
verbatim as `LICENSE` and `upstream/LICENSE`.

References:

- Jorjani, H. et al. (2003), *A Simple Method for Weighted Bending of Genetic
  (Co)variance Matrices*, Journal of Dairy Science 86:677-679.
- Schaeffer, L. R. (2014), *Making covariance matrices positive definite*.
