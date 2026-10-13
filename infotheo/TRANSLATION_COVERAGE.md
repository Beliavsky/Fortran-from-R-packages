# Computational API coverage

Status: **substantial**. Basis: **exported computational R functions**.
Mapped: **8 of 8 (100%)**. Unmapped: none. A mapping count is not an assertion
of complete R-interface compatibility. NAMESPACE exports names by pattern;
there are exactly eight exported R functions in the two upstream R files.

| R function | R source | Fortran module | Public procedure(s) | Status |
|---|---|---|---|---|
| entropy | upstream/R/entropy.R | infotheo | entropy | substantial |
| condentropy | upstream/R/entropy.R | infotheo | condentropy | substantial |
| mutinformation | upstream/R/entropy.R | infotheo | mutinformation, mutinformation_matrix | substantial |
| condinformation | upstream/R/entropy.R | infotheo | condinformation | substantial |
| multiinformation | upstream/R/entropy.R | infotheo | multiinformation | substantial |
| interinformation | upstream/R/entropy.R | infotheo | interinformation | substantial |
| discretize | upstream/R/discretize.R | infotheo | discretize | substantial |
| natstobits | upstream/R/entropy.R | infotheo | natstobits | complete |

All procedures are in `src/infotheo.f90`. `entropy` is a generic for integer
rank-one/rank-two data; remaining information functions accept rank-two blocks.
The four upstream estimators and three discretization methods are implemented.
The pairwise matrix operation is included, not just scalar mutual information.
Missing-data support is retained via an explicit integer sentinel.

Material differences: corrected SG prior with `sg_r` compatibility option;
correct full-block conditional information; full-precision nats conversion;
finite/NA/extreme-range safety; singleton shrink entropy zero; an exponential
interaction guard of 20 variables. Missing-value deletions remain per entropy
term and scalar versus matrix MI truncation follows R. Fortran arrays replace
R formula/data-frame/factor coercion and naming; no R object construction or
R warning side effects are supplied. See README.md for exact contracts and
NOTICE.md for provenance and dependency reuse.
