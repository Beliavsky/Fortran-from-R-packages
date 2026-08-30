# Translation notes

## Scope

The translation targets the non-R computational behavior of `proxy` 0.4-29:

- built-in distance and similarity formulas from `R/dissimilarities.R` and `R/similarities.R`;
- optimized missing-data semantics from `src/distance.c` for Euclidean, Manhattan, Minkowski, Canberra, supremum, cosine, extended Jaccard/Dice, and fuzzy Jaccard;
- matrix auto/cross/pairwise application loops from `src/apply.c`;
- packed `dist` indexing, subsetting, row sums/means, and row/column index helpers from `src/util.c`;
- mixed-variable Gower preprocessing and scoring;
- Levenshtein edit distance, implemented natively rather than through the optional R `cba` dependency;
- an extensible native numeric/binary registry replacing the R environment/S3 registry machinery.

R formula parsing, S3 printing/class attributes, data-frame/list dispatch, character subscript emulation, and namespace machinery are language-interface concerns and are not translated literally.

## Native API decisions

Fortran numeric matrices are observation-by-variable. Auto/cross/pairwise modes are separate explicit routines instead of being selected by nullable R arguments. Callers that conceptually want upstream `by_rows = FALSE` can transpose the input before calling the native API.

Mixed Gower variables use numeric encodings plus explicit type codes. This preserves the statistical computation while avoiding an R-like dynamically typed data frame in the Fortran core. Unused R factor levels are metadata and therefore are not represented by the numeric factor-code API.

Mahalanobis auto-distance estimates `cov(x)` when no covariance is supplied. Cross Mahalanobis follows upstream `cov(x,y)` behavior when `x` and `y` have equal row counts; otherwise callers must supply a covariance matrix.

Levenshtein is implemented natively. Its scalar and character-array routines are separate from the real-matrix dispatcher because Fortran is statically typed.

The current typed custom registry supports numeric-vector callbacks and binary contingency-count callbacks. It does not attempt to reproduce the arbitrary heterogeneous R-function registry or its S3 metadata/validation surface.

## Literal upstream quirks retained

Several behaviors are deliberately translated as the package implements them rather than replaced by more conventional textbook formulas:

- `WaveHedges` uses R scalar `min(x,y)` and `max(x,y)` semantics.
- `Soergel` divides by the scalar global `max(x,y)` used in upstream R code.
- exact infinite Minkowski exponents retain the behavior of the upstream C fast path.
- the matrix fast path for extended Dice uses `2*xy/(xx+yy)`, matching upstream C even though the R fallback function omits the factor two and is marked `FIXME` upstream.
- binary Jaccard returns unit similarity for two all-false vectors on the matrix fast path.
- formulas with zero denominators preserve upstream NaN/Inf behavior without relying on trapped hardware division by zero.

## Validation

The strict suite exercises every registered built-in family, auto/cross/pairwise consistency, missing-value compensation, all binary coefficients, nominal association measures, mixed Gower data, custom registry callbacks, Levenshtein distance, and packed `dist` helpers.

Additional fixtures are copied from upstream saved test output (`tests/apply.Rout.save` and `tests/util.Rout.save`) for cosine auto/cross similarities, packed row sums, and `row.dist`/`col.dist` indices.

No R executable is installed in the translation environment, so live R-vs-Fortran execution was not possible. The retained upstream snapshot is included for source-level parity review.

## Source conventions

- `dp = real64` is defined exactly once in `proxy_kinds` and re-exported by `proxy`.
- Every maintained dummy argument has explicit `INTENT` or `VALUE`.
- Every dummy argument is declared on its own line.
- Every dummy declaration has a meaningful trailing FORD `!!` documentation comment.
- Maintained source is free-form, stays within the normal 132-column source limit, and is formatted to be compatible with `fprettify`.
- `fprettify` was not installed in the release-validation environment, so compatibility was enforced by source style/audit rather than by running that executable.
