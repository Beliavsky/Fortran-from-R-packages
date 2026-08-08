# Translation coverage and source differences

## Directly translated computational code

### `src/mcga.c` and exported R `mcga()`

Translated into `mcga_engine.f90` and `mcga_bytes.f90`:

- random population generation;
- scalar cost handling;
- elitist sorting;
- pairwise tournament selection;
- uniform byte crossover;
- byte mutation;
- population replacement;
- final cost sorting.

The exported R wrapper does not evaluate the randomized initial population
before the first tournament. The Fortran `mcga_optimize` intentionally matches
that behavior: initial costs are zero and the first tournament is therefore
effectively random apart from tie behavior.

The bounds `minval,maxval` are only initialization bounds in upstream `mcga()`;
they are not enforced after raw byte operators. The Fortran engine preserves
that behavior.

### `src/multi_mcga.c` and exported R `multi_mcga()`

Translated including the original rank-score algorithm. For chromosome `i`,
the score is incremented once for each chromosome `j` for which `i` has a
strictly smaller value in at least one objective. This is the actual upstream
algorithm, although the documentation calls it non-dominated sorting.

The upstream tournament code uses `>` for the first parent rank comparison and
`<` for the second. The Fortran native multi-objective engine preserves this
asymmetry exactly.

### `src/typeconversations.cpp`

Translated:

- `MaxDouble`, `SizeOfDouble`, `SizeOfInt`, `SizeOfLong`;
- double/byte conversions;
- one-point, two-point, and uniform crossover;
- +/-1 byte mutation with 0/255 wraparound;
- random byte mutation;
- bounds repair.

`TRANSFER` replaces C++ pointer reinterpretation. This is the correct Fortran
mechanism for retaining the host bit pattern. Consequently, like upstream,
results may be machine/endianness dependent.

### `R/oplibrary.R`

All mcga-specific mutation and crossover operators are translated.

Three R functions contain apparent implementation/indexing defects. The
Fortran routines implement their documented/intended operators rather than
reproducing malformed R behavior:

1. `sbx_crossover`: upstream uses the final scalar `betaq[i]` in part of each
   vector expression; Fortran uses the standard elementwise SBX formula.
2. `linear_crossover`: upstream uses `[[setdiff(...)]]`, which is recursive R
   indexing rather than selecting the remaining two list elements; Fortran
   evaluates all three candidates and returns the two highest-fitness ones.
3. `unfair_average_crossover`: upstream loops only through `1:j`, making its
   `else` branch unreachable and leaving later genes zero; Fortran loops over
   all genes and applies the two documented branches around `j`.

These are deliberate correctness fixes, not accidental translation drift.

## External dependency not reimplemented

`mcga2()` contains no independent optimizer: it calls `GA::ga` from the
external `GA` package while supplying mcga's byte operators. Since the GA
package was not part of the supplied source archive, its population generator,
selection method, stopping rules, parallel evaluation, S4 result type, and
monitoring system are outside this translation. The complete mcga operator
library needed by such a driver is included.

## Omitted non-computational code

- R `.Call`/Rcpp wrappers and registration;
- R environments and SEXP handling;
- printing/dump helpers whose purpose is diagnostic display;
- Rd/vignette plotting examples;
- GA S4 object plumbing and monitoring.

The original package is retained under `original/mcga-master/`.
