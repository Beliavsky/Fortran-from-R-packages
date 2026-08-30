# Translation notes

## Upstream baseline

- Package: `rpart`
- Version: 4.1.27
- Upstream package date: 2026-03-26
- License: GPL-2 | GPL-3
- Authors listed by upstream: Terry Therneau, Beth Atkinson, Brian Ripley

The complete source package used for this translation is retained in
`upstream/`.

## Source mapping

| Upstream computational area | Fortran implementation |
| --- | --- |
| `anova.c` | `rpart_methods.f90` |
| `gini.c`, `graycode.c` | `rpart_methods.f90` |
| `poisson.c` | `rpart_methods.f90` |
| `partition.c`, `bsplit.c` | `rpart_tree.f90` |
| `surrogate.c`, `choose_surg.c`, `nodesplit.c` | `rpart_tree.f90` |
| `make_cp_list.c`, `make_cp_table.c`, `fix_cp.c` | `rpart_cp.f90`, `rpart_tree.f90` |
| `xval.c` | `rpart_fit.f90` |
| `xpred.c`, `xpred.rpart.R` | `rpart_xpred.f90` |
| `pred_rpart.c`, `predict.rpart.R` | `rpart_predict.f90` |
| `rpart.exp.R`, `rpartexp2.c` | `rpart_survival.f90` |
| `rpart.control.R` | `rpart_control_api.f90`, `rpart_utils.f90` |
| `importance.R` | `rpart_cp.f90` |
| `prune.rpart.R` | `rpart_cp.f90` |

## Design choices

The C implementation relies on package-global state, R's memory manager, and
function-pointer tables. The Fortran translation stores method state in
`rpart_model` and passes it explicitly through typed procedures. Recursive
nodes are represented by allocatable derived-type children.

The original code's categorical Gray-code search is retained. Binary
classification and regression/Poisson categorical predictors use the same
ordered-category reductions as upstream; multiclass unordered classification
uses full Gray-code subset enumeration.

Missing predictors are represented by any non-finite `real(dp)` value and are
routed by primary/surrogate/majority logic corresponding to rpart.

## CP and cross-validation

Node complexities are stored internally in risk units. Public CP-table values
are relative to the root risk, as in R's printed `cptable`.

Fold trees rescale the full-tree CP by the training-weight fraction and the
fold root risk, matching the scaling in upstream `xval.c`/`xpred.c`.
Classification fold trees retain full-data prior/frequency scaling, while
Poisson/exponential fold trees recompute their empirical-Bayes rate prior, as
the upstream initializers do.

The R RNG is not reproduced. Pass explicit fold groups when fold identity
needs to match another implementation.

## Survival preprocessing

`rpart_exp_transform_right` and `rpart_exp_transform_startstop` reproduce the
current R `rpart.exp` preprocessing: event-time interval construction,
roundoff-scale event-time amalgamation, the 1000-interval cap, cumulative
hazard stretching, and optional offset multiplication. The current upstream R
implementation passes case weights to its helper but does not use them in the
interval-rate calculation; the Fortran translation therefore does not apply
weights during this transformation either.

## Deliberately omitted R-only surfaces

These are not missing tree algorithms:

- formula/model-frame construction and factor names/labels;
- S3 print/summary/text methods;
- `plot.rpart`, `plotcp`, `post`, `path.rpart`, `snip.rpart`, and graphics;
- R's NA action object bookkeeping;
- the R callback bridge for `method="user"`.

The last item is an extension mechanism rather than one of rpart's built-in
statistical methods. A future native Fortran callback API could provide an
analog without carrying R runtime machinery into this package.

`meanvar.rpart`, `rsq.rpart`, and `roc.rpart` are primarily reporting/plotting
helpers and are not dedicated APIs in v0.1.0; their required fitted values,
CP-table values, node responses, and class probabilities are available from
the model and prediction APIs.

## Native-interface qualifications

- Direct mutation of an `rpart_control` object cannot emulate R's
  missing-argument logic for coupled `minsplit`/`minbucket` defaults. Use
  `rpart_make_control` when that behavior matters.
- `prune_model` cannot recompute training `where`/`fitted` without retaining
  the original design matrix, so it deallocates those two cached arrays after
  structural pruning. Predictions from the pruned tree remain available via
  the prediction API.
- Ordered R factors should be passed as their ordered numeric codes with
  `ncat(j)=0`, which is how upstream sends ordered factors to its C kernel.

## Precision and source policy

`rpart_kinds.f90` defines `dp = real64` exactly once. Maintained source imports
that same kind throughout and uses `_dp` constants. No maintained Fortran uses
`double precision`, `real*8`, `kind(0.0d0)`, or D-exponent literals.
