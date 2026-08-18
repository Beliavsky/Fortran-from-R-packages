# Translation notes

## Scope

This is a computational translation of `hyper2` 3.2-3 to modern Fortran. The
R package represents likelihoods through Rcpp/STL maps plus R S3/S4 classes;
the Fortran translation replaces those with explicit derived types and
canonical sparse arrays.

### `hyper2` likelihood algebra

`hyper2_model` stores each bracket as a sorted set of integer player IDs and a
real power. Equal brackets are combined at insertion time. The implementation
covers the work performed upstream by `identityL`, `addL`, access/overwrite,
assignment, equality, evaluation, differentiation, and low-level Hessian code.

`loglik`, `gradient`, `gradient_full`, and `hessian_independent` accept full
simplex vectors. Gradients/Hessians returned by the short forms are with respect
to the first `n-1` probabilities with the last probability supplied by the
simplex constraint.

### `hyper3`

`hyper3_model` stores a nonnegative weight for each member of a bracket. This
replaces the upstream map-of-maps Rcpp representation. Weighted evaluation,
gradients, Hessians, conversions, weight changes, home/away/draw constructors,
and repeated-player race likelihoods are native Fortran. `hyper3_matrix`
provides the matrix constructor corresponding to upstream `hyper3_m`; callers
can also construct arbitrary weighted brackets directly with `%add_term`.

### Maximum likelihood

The upstream package delegates constrained maximization to R optimizers
(`constrOptim`, `alabama`). `maxp` uses exact gradients/Hessians with safeguarded
Newton/backtracking steps on the simplex, and supports linear inequality
constraints. `maxp_simplex` provides multistart fitting. `maxp_lsl` uses a
softmax parameterization and numerical derivatives for list-of-sum likelihoods.

### Grouped ranks and list likelihoods

`general_grouped_rank_likelihood` expands tied groups with the translated
`partitions::perms` algorithm and returns a `hyper2_list`. All alternatives are
remapped onto one canonical player order before evaluation. `suplist_add`
implements the Cartesian product semantics of upstream `Ops.suplist`, and
`suplist_scale` implements repeated-product integer powers rather than merely
scaling term exponents.

### Integration

Upstream `B()` calls `cubature::adaptIntegrate`. `hyper2_B` uses the vendored
pure-Fortran `cubature-fortran` translation with the same stick-breaking
simplex transform and Jacobian. `dhyper2`, `mgf`, and `mean_hyper2` are built on
that normalizer. The callback context is module-saved and therefore a single
`hyper2_B` call is not re-entrant/thread-safe.

### Simulation

The R package's `rdirichlet`, rank/race generators, and support Metropolis
sampler are translated natively. Gamma variates use Marsaglia-Tsang sampling.
Fortran's intrinsic RNG is used, with `seed_rng()` providing deterministic
seeding within a compiler/runtime implementation.

## Computational routines represented directly

- hyper2/hyper3 construction and algebra
- log-likelihood/evaluation and derivatives
- exact simplex Hessians
- maxp/maxp_simplex
- dirichlet/GD/GD_wong/dirichlet3
- rankvec_likelihood/ordervec2supp/ordervec2supp3/ordervec2supp3a
- pairwise, home_away, home_away3, home_draw_away3
- Zermelo iteration
- as_hyper3/hyper3_to_hyper2, setweight, pwa/pwa3/pwa23
- keep/discard/substitute/balance
- grouped-rank `ggrl` equivalent, suplist, lsl
- B/dhyper2/mgf/mean_hyper2
- rank-table/order-table conversions and support
- rdirichlet/rp/rank/race/hyper2/hyper3 simulation primitives
- equalp.test and knownp.test likelihood-ratio calculations

## R-specific or presentation infrastructure omitted

The following are not meaningful as direct Fortran APIs and are omitted:

- S3/S4 class registration, replacement methods, magrittr syntax, disord/frab
  hashing infrastructure, and Rcpp registration;
- printing/summary formatting and wiki/data-frame coercion helpers;
- `ordertransplot` and other plotting-only behavior;
- package datasets and documentation-only convenience objects;
- CRAN parallel/process plumbing.

A few highly R-interface-specific convenience wrappers (for example formula-1
points-system list construction, attempt-table data-frame formatting, and the
one-sided/profile variants of the hypothesis-test wrappers) are not reproduced
as one-for-one APIs. Their underlying likelihood construction, fitting,
simulation, and LR-test machinery is available through the translated modules.

## Fortran-specific correctness fixes

Fortran does not guarantee short-circuit evaluation of `.and.`. Runtime bounds
checking exposed two places where a C/R-style short-circuit assumption would
have allowed an index-zero access during weighted bracket canonicalization and
array comparison. They are written as nested tests in this translation.

Grouped-rank alternatives may enumerate players in different order. The R
objects carry names, whereas a numeric Fortran probability vector is positional.
The translation therefore remaps every alternative to a common canonical player
order before a suplist is evaluated.

## Source form

All compiled Fortran is free-format `.f90`, ASCII, and intended to stay within
the standard 132-column free-form limit. The FPM manifest disables implicit
typing and implicit external interfaces.
