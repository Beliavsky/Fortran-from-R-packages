# infotheo: information-theoretic measures in modern Fortran

A substantial translation of Patrick E. Meyer's R **infotheo 1.2.0.1**.
It provides joint/conditional entropy, mutual and conditional information,
pairwise mutual-information matrices, total correlation, co-information,
three discretizers, four entropy estimators, and conversion from nats to bits.
All maintained Fortran is free form. Every library procedure is pure;
`natstobits` is also elemental. GPL-3.0-or-later; see `LICENSE` and `NOTICE.md`.

## Build and run

Place this directory at the root of
[Beliavsky/Fortran-from-R-packages](https://github.com/Beliavsky/Fortran-from-R-packages),
alongside the existing `rfortran-core` directory, then run:

```sh
cd infotheo
fpm build
fpm test
fpm run --example
```

The sole dependency is `rfortran-core = { path = "../rfortran-core" }`.
Its existing `dp`, `r_digamma`, stable `r_order`, and missing-aware `r_sort`
are reused. No dependency implementation, BLAS, LAPACK, R runtime, C++ compiler,
or external numerical library is bundled or needed. The delivered R/C++ sources
under `upstream/` are provenance only and are not compiled by FPM. The Fortran
build uses portable FPM/gfortran facilities, including on Windows; native Windows
execution has not been performed. GNU Fortran 14.2.0/FPM 0.13.0 Linux validation
is recorded in `VALIDATION.md`.

Do not use `-ffast-math`, `-Ofast`, or `-ffinite-math-only`: missing-data behavior
and explicit NaN error results require IEEE semantics.

## API and data conventions

Import the public module `infotheo`; it re-exports `dp` from `r_kinds`.
Discrete data use default-integer category labels, observations in rows and
variables in columns. Any integer label is allowed except `na_integer`, the
reserved missing sentinel `-huge(0)-1`. The tested gfortran default integer kind
is 32 bits, and the shared `dp` is binary64. Factors/strings must be encoded
before calling these routines. R coercion, formula/data-frame handling,
column names, warnings, and object construction are not translated.

```fortran
use infotheo, only : dp, entropy, condinformation, mutinformation_matrix, natstobits
integer :: x(4,3)
real(dp) :: h
real(dp), allocatable :: mi(:,:)
x(:,1) = [0, 0, 1, 1]
x(:,2) = [0, 1, 0, 1]
x(:,3) = [0, 1, 1, 0]
h = natstobits(entropy(x))  ! 2 bits
h = natstobits(condinformation(x(:,1:1), x(:,2:2), x(:,3:3)))  ! 1 bit
mi = mutinformation_matrix(x)
```

| Procedure | Result and inputs |
|---|---|
| `entropy(x, method)` | Scalar joint entropy; `x` may be an integer vector or matrix |
| `condentropy(x, y, method)` | H(X,Y) - H(Y); integer matrices; optional `y` defaults to H(X) |
| `mutinformation(x, y, method)` | max(0, H(X) + H(Y) - H(X,Y)); integer matrix blocks |
| `mutinformation_matrix(x, method)` | Allocatable square matrix of pairwise MI; diagonal H(X_i) |
| `condinformation(x, y, s, method)` | H(S,X) + H(S,Y) - H(S) - H(S,X,Y); optional `s` defaults to scalar MI |
| `multiinformation(x, method)` | Sum of marginal entropies minus joint entropy |
| `interinformation(x, method)` | Alternating nonempty-subset entropy sum, starting with positive singleton terms |
| `discretize(x, bins, disc, nbins, status)` | Pure subroutine; real(dp) matrix to allocatable integer matrix of the same shape |
| `natstobits(h)` | Elemental H/log(2), preserving shape for real(dp) arrays |

Except `entropy`, discrete inputs are rank-two matrices. To pass one column,
use `x(:,j:j)` or `reshape(vector,[size(vector),1])`.
All information values are in **nats** before conversion. Optional estimator
`method` is case-sensitive, defaults to `"emp"`, and accepts:

- `"emp"`: Shannon entropy of observed frequencies
- `"mm"`: empirical entropy plus (K - 1)/(2N), using K observed states
- `"sg"`: posterior expected Dirichlet entropy with symmetric prior beta = 1/K
- `"shrink"`: the upstream shrinkage Dirichlet estimator, including the
  uniform-target limit. This is not the plug-in entropy of shrunk frequencies.
- `"sg_r"`: compatibility extension reproducing upstream's integer-division
  prior, beta = 0 for K > 1. Use this to compare with R's `method="sg"`.

For posterior parameters alpha_i = count_i + beta, the Dirichlet entropy is
sum_i alpha_i / sum(alpha) * [digamma(sum(alpha)+1) - digamma(alpha_i+1)].
For shrinkage, write w = sum_i (count_i/N)^2 and use
lambda = K(1-w)/[(N-1)(Kw-1)], bounded to [0,1]. If lambda < 1,
beta = lambda N/[(1-lambda)K]; otherwise entropy is log(K).
Only observed states are included, as in upstream.

### Missing data, invalid arguments, and computational limits

Each entropy removes observations missing in any of its participating columns.
Combined quantities retain R's **separate available-case calculation** for each
entropy term. With missingness, those terms can use different rows; negative
conditional information and failures of usual nonnegativity identities can
result. Filter to common complete rows yourself when that is the desired
statistical convention. Scalar MI clips negative results to zero, while matrix
MI, conditional quantities, total correlation, and interaction values do not,
matching the distinct R operations.

No complete observations, unknown estimators, row-count mismatches, or an
interaction request above `max_interaction_variables = 20` return a quiet NaN.
A nonempty observation set with zero variables has entropy zero. Constant or
single-observation samples have entropy zero for every estimator, including
shrinkage. Interaction needs 2^p-1 entropy evaluations; the 20-column guard
replaces the unrealistic upstream 500-column guard. Counting uses shared stable
sorting, O(p N log N) per entropy with O(N) auxiliary storage, not a quadratic
search for distinct rows. Allocation failure uses normal Fortran runtime behavior.

### Discretization

`call discretize(x, bins, disc, nbins, status)` accepts:

- `disc="equalfreq"` (default): per-column sorted cutpoints, preserving the
  upstream rank calculation and assigning tied values to the earliest matching
  bin. For N=8 and B=3 the internal cut ranks are 3 and 5. The historical
  remainder rule is retained even when bin counts are not optimally balanced.
- `disc="equalwidth"`: per-column, left-closed intervals; the observed maximum
  belongs to the final bin. Constant columns map to bin 1.
- `disc="globalequalwidth"`: a common range over all columns, right-closed
  intervals as in R `cut`; constant data map to ceiling(B/2).

`nbins` defaults to `int(real(nrows,dp)**(1.0_dp/3.0_dp))`, matching R's
floating exponent followed by integer conversion. Explicit counts must be
positive; global width requires at least two. Equal frequency also requires
`nbins <= number of nonmissing observations` for each nonempty column. A
missing-only column is allowed and remains missing. NaNs become `na_integer`;
infinities are invalid. Output labels are 1:nbins and may leave empty bins.
`status` is optional and is `infotheo_ok=0` on success, `infotheo_invalid=1`
on invalid arguments. On failure, the output is allocated with the input shape
and entirely filled with `na_integer`. An empty-row input has an invalid default
bin count; an explicit valid count can produce an empty output successfully.

## Deliberate corrections and compatibility boundaries

1. **SG prior:** upstream `1/freq.size()` performs integer division, changing
   the documented estimator. `sg` uses floating 1/K; `sg_r` retains its numerical
   behavior for valid samples. Singleton shrinkage is defined as zero instead
   of upstream's 0/0 NaN; square calculations use real probabilities, avoiding
   upstream 32-bit overflow for large samples.
2. **Multicolumn conditioning:** upstream `condinformation` selects hardcoded
   columns 1:2 and columns 1,3. This fails for random-vector blocks. The Fortran
   API uses all X, Y, and S columns according to the documented identity.
3. **Nats conversion:** exact `1/log(2)` replaces the truncated `1.442695`.
4. **Discretizer safety:** real extrema replace upstream integer extrema;
   missing markers remain exactly `na_integer` instead of being incremented;
   sorting explicitly removes NaNs; invalid counts are rejected before unsafe
   indexing; huge finite ranges avoid subtraction overflow. Missing-only global
   input returns missing output rather than failing inside R `cut`. This does
   not reproduce unsafe source behavior or numerical accidents at sub-ULP cuts.
5. **R interfaces:** the matrix MI operation is a separate named function.
   This also avoids the R wrapper's `is.atomic(matrix)` rejection: upstream
   reference tests pass a data frame to that operation.

See `VALIDATION.md` for independently checked formulas, actual-R fixtures, and
edge cases. `tools/generate_reference.R` regenerates the checked-in actual-R
fixtures when R infotheo 1.2.0.1 is installed; R is not needed for ordinary tests.

## Translation coverage

Package status: **substantial**. **8 of 8 (100%)** exported computational R
functions are mapped. This fraction measures mapped computational functions,
not complete R compatibility. Internal C/C++ entry points and R-only interfaces
are not separate functions in the denominator. There are no unmapped exported
computational functions.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `entropy` | `infotheo` | `entropy` | substantial |
| `condentropy` | `infotheo` | `condentropy` | substantial |
| `mutinformation` | `infotheo` | `mutinformation`, `mutinformation_matrix` | substantial |
| `condinformation` | `infotheo` | `condinformation` | substantial |
| `multiinformation` | `infotheo` | `multiinformation` | substantial |
| `interinformation` | `infotheo` | `interinformation` | substantial |
| `discretize` | `infotheo` | `discretize` | substantial |
| `natstobits` | `infotheo` | `natstobits` | complete |

The machine-readable mappings in `fpm.toml` and `TRANSLATION_COVERAGE.md`
use the same basis and counts. The original requested translation instructions
are preserved verbatim in `TRANSLATION_PROMPT.txt`.
