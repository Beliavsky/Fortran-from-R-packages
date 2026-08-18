# Translation notes

## Scope

This release translates the computational code of `partitions` 1.10-9 to
modern Fortran 2018.  All Fortran source is free format (`.f90`).

Translated functionality includes:

- unrestricted additive partitions and successor/first/last operations;
- partitions into distinct parts;
- partitions into a fixed number of parts, with or without zero padding;
- bounded/block partitions, including the upstream enumeration order;
- unrestricted and fixed-length compositions;
- partition-count functions corresponding to R `P()`, `Q()`, `R()`, and
  `S()`;
- conjugate partitions and Durfee-square size;
- set partitions with specified block sizes;
- restricted set-partition ordering matrices;
- lexicographic permutations and Knuth plain-change permutations;
- unique multiset permutations and selections;
- multinomial permutations, all-binomial selections, generalized riffles,
  and two-packet riffles;
- integer/binary and composition/binary conversion.

## API organization

The main module is:

```fortran
use partitions
```

The translation provides two naming layers.

The explicit modern API uses names such as:

- `partition_count`, `distinct_partition_count`
- `restricted_partition_count`, `block_partition_count`
- `distinct_parts`, `restricted_parts`, `block_parts`
- `first_part`, `next_part`, `is_last_part`
- `set_partitions`, `multiset_permutations`
- `multinomial_permutations`, `generalized_riffles`
- `to_binary`, `to_decimal`, `composition_to_binary`

A compatibility layer also exposes familiar package names where practical:

- `p`, `q`, `r`, `s`
- `parts`, `diffparts`, `restrictedparts`, `blockparts`, `compositions`
- `firstpart`, `nextpart`, `islastpart`, and corresponding variants
- `setparts`, `restrictedsetparts`, `restrictedsetparts2`
- `perms`, `plainperms`, `mset`, `multiset`, `multinomial`
- `allbinom`, `genrif`, `riffle`
- `tobin`, `todec`, `comptobin`, `bintocomp`

R's `P(..., give=TRUE)` and `Q(..., give=TRUE)` rank-changing behavior is
represented by the separate functions `partition_numbers()` and
`distinct_partition_numbers()`.

## Representation differences

R represents enumerations as matrices whose columns are partitions.  The
Fortran translation deliberately retains that convention: results are
allocatable integer matrices with one object per column.

R class, print, summary, row-name, and `sets`/`equivalence` object machinery
is not reproduced.  In particular, `print.partition`, `summary.partition`,
`print.equivalence`, `listParts`, `vec_to_set`, `vec_to_eq`, and `condense`
are presentation/object-layer code rather than numerical enumeration code.
The underlying set-partition membership labels and restricted ordering
matrices are available directly.

## Numerical changes

The upstream C routines return partition counts through double precision and
use 64-bit integer work arrays.  The Fortran translation returns exact
`integer(int64)` counts while they fit.  Overflow is detected and reported
instead of silently wrapping.

The R implementation uses `gmp` factorials to calculate the number of set
partitions before allocating the output.  The Fortran implementation cancels
factorials by prime exponents, avoiding unnecessary intermediate factorial
overflow without adding a GMP dependency.  The final enumeration size must
still fit `int64` and available memory.

The upstream `S()` uses the `polynom` package.  The Fortran implementation
computes the same coefficient by exact integer dynamic programming.

## Algorithm fidelity

The native algorithms were translated directly where ordering is observable:

- `c_nextpart`
- `c_nextdiffpart`
- Andrews' restricted-partition successor
- `c_nextblockpart`
- Knuth Algorithm L lexicographic permutations
- Knuth Algorithm P plain changes
- the C++ `Partitions` recursion used by `setparts()`

The set-partition translation reproduces the exact upstream column order for
`setparts(c(2,1,1))`, including the regression case in the R test suite.

## Validation

The included tests cover the upstream fixed reference values

- `P(100) = 190569292`
- `Q(100) = 444793`
- `R(5,12) = 13`
- `S(rep(1:4, each=2), 5) = 474`

as well as enumeration/successor agreement, partition sums and ordering,
conjugation/Durfee examples, exact set-partition order, compositions,
binary conversion, permutations, multisets, multinomial arrangements, and
riffles.
