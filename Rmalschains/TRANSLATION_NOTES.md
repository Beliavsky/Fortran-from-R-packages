# Translation notes

## Source

Translated from the attached CRAN source tree:

* Package: Rmalschains
* Version: 0.2-11
* Title: Continuous Optimization using Memetic Algorithms with Local Search
  Chains (MA-LS-Chains)
* Declared package license: GPL-3

The R package adapts Daniel Molina's librealea implementation.

## R-to-Fortran mapping

`malschains.control()` maps to `type(mals_control)` and the convenience function
`malschains_control()`.

`malschains()` maps to `malschains_optimize()`.

R closures/environments are replaced by an explicit Fortran procedure callback:

```fortran
function objective(x) result(value)
  real(dp), intent(in) :: x(:)
  real(dp) :: value
end function
```

The returned R list maps to `type(mals_result)`.

## SSGA

The translated evolutionary phase follows the exact active wrapper setup:

* `CrossBLX(alpha)`
* `MutationBGA()` wrapped in `Mutation` with default probability 0.125
* `SelectNAM(3)`
* `ReplaceWorst()`

`SelectNAM(3)` chooses a random first parent and selects the most distant of
three separately sampled candidates as the second parent.

The BGA mutation uses the 16-term binary geometric increment from
`MutationBGA::mutate()` and a perturbation range of 10% of the variable domain.

## MA-LS-Chains scheduling

The C++ routine

```text
calculateFrec(nevalalg, nevalls, intensity, effort)
```

is translated directly.  This is important: `effort` is enforced against
**cumulative** EA/LS evaluation totals, not by applying a fixed phase fraction
at every cycle.

Per-individual LS state is retained only when the local search improves the
individual by the package's `hasImprovedEnough()` criterion.  Failure marks the
individual as `non_improved`, discards its LS state, and removes it from the
normal LS-candidate set until it is later replaced by SSGA.

The wrapper disables restart (`ma->setRestart(NULL)`), so the Fortran R-facing
path does not invent a restart operation.

## Local searches

### Solis-Wets (`sw`)

Translated from `solis.cc`, including bias adaptation, paired forward/reverse
evaluations, success/failure counters, and multiplicative delta adjustment.

### Sparse Solis-Wets (`ssw`)

Translated from `solisn2.cc` strategy 3: each LS application marks dimensions
with probability 0.1, and only those dimensions receive Gaussian perturbations
and delta adaptation.

The original initialization reads an uninitialized `delta` element while
clamping `delta_init`; the port uses the intended `delta_init` value instead.

### Simplex

Translated from `simplex.cc`: the initial simplex consists of the current point
plus one point shifted by 10% of the range in every dimension.  Reflection,
expansion, contraction and the original contraction-triggered simplex shrink
are retained.  The simplex is persistent across chain calls.

### MTS1/MTS2

The native selector accepts these methods even though the R help emphasizes
four local searches.  Their coordinate/random-subset steps and persistent delta
state are translated.

### CMA-ES

The local CMA-ES chain retains persistent mean, covariance, eigenbasis,
evolution paths and step size.  Initial coordinate standard deviations follow
the package's neighborhood idea: half the nearest nonzero coordinate distance
plus `0.001`.

The historical package embeds Nikolaus Hansen's large C CMA-ES implementation.
The Fortran port implements the same standard covariance/path/step-size
adaptation equations directly and uses LAPACK `DSYEV`; it does not preserve the
exact internal C data layout or random stream.

`lsParam1` and `lsParam2` map to CMA-ES `lambda` and `mu`.

## Evaluation accounting

Rmalschains has multiple counting layers.  In particular, initial-population
values supplied by R are evaluated after MA initialization with `Problem::eval`
and therefore bypass the `Running` counter.

`mals_result` exposes the package-style counts:

* `num_eval_ea`
* `num_eval_ls`

and a new diagnostic:

* `actual_nfe`

for every Fortran callback invocation.

## R infrastructure intentionally omitted

* Rcpp `.Call` registration and SEXP conversion;
* R environments and compiled-R callback handling;
* S3 `print.malschains` formatting;
* console debug macros and R printing;
* serialization/distributed-computing hooks from dormant librealea classes;
* inactive vendor algorithms not reachable from `malschains()`.
