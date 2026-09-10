# API coverage

The coverage basis is the 13 exported computational functions in the upstream `NAMESPACE`, which exports all names. Nested helpers are not counted separately.

| R function | Computational coverage | Material differences |
| --- | --- | --- |
| `combn` | Combination generation and ordering | No `fun` callback, R list/simplify machinery, character arrays, or scalar-to-`seq` coercion. |
| `combn2` | Pair-index and numeric pair generation | Explicit Fortran overloads replace R missing-argument dispatch and generic vectors. |
| `dmnom` | Multinomial mass, recycling, normalization, default `size` | Numeric arrays only; no R warnings/coercions. |
| `fact` | `gamma(x+1)` | Elemental real Fortran interface. |
| `hcube` | Lattice ordering, scale, translation, recycling | Invalid dimensions return status rather than R errors. |
| `logfact` | `lgamma(x+1)` | Elemental real Fortran interface. |
| `nCm` | General real `n`, integer-valued `m`, reflection, product branch, gamma branch, rounding, recycling | Impractically large `m` outside default-integer loop bounds returns NaN. |
| `nsimplex` | Point counts and recycling | Numeric Fortran arrays instead of R vectors. |
| `permn` | Minimal-change permutation algorithm and order | No callback/list interface or arbitrary R types. |
| `rmultinomial` | Varying `n`, varying probability rows, recycling | Explicit Park-Miller RNG state; not R-stream compatible. |
| `rmultz2` | Varying `n`, fixed probabilities, upstream output orientation | Explicit Park-Miller RNG state; not R-stream compatible. |
| `x2u` | Bin-index expansion and integer/real labels | Counts are integer and labels are integer or real rather than arbitrary R objects. |
| `xsimplex` | Simplex point generation and order | No callback or R list/simplify interface. |

The package-level status is **substantial** even though all 13 computational functions are mapped, because R's dynamic callback/object interfaces and RNG stream are intentionally not reproduced.
