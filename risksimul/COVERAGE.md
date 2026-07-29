# Computational coverage

Upstream package: `riskSimul` 0.1.2.

## Translated routines

| Upstream routine | Fortran implementation | Status |
|---|---|---|
| `OrthMat` | `orthogonal_completion` | Implemented with stable Gram-Schmidt completion |
| `OptAllocHeur` | `optimal_allocation_heuristic` | Implemented; corrected and compatibility modes |
| `bvec` | internal multiresponse relative-error weights | Implemented |
| `amat` | internal stratum-variance matrix construction | Implemented |
| `ReturnCopula` | `return_copula`, `portfolio_return_one` | Implemented |
| `TailLossProb` | `tail_loss_response` | Implemented |
| `Excess` | `excess_response` | Implemented |
| `Touch` | `touch_value` | Implemented |
| `Alg2` | `algorithm_2` | Implemented |
| `Alg3` | `algorithm_3` | Implemented with native pattern search |
| `NVCopula` | `naive_copula` / `NVTCopula` | Implemented |
| `searchx` | `search_threshold` | Implemented |
| `NVCopulaMT` | vector-threshold `naive_copula` | Implemented |
| `SISCopula` | scalar-threshold `stratified_copula` | Implemented |
| `SISCopulaMT` | vector-threshold `stratified_copula` | Implemented |
| `NVTCopula` | `NVTCopula` | Implemented |
| `SISTCopula` | `SISTCopula` | Implemented |
| `new.portfobj` | `new_portfolio` / `new_portfobj` | Implemented |

## Marginal models

- Student-t marginals: complete.
- Generalized-hyperbolic marginals: complete numerical support through native
  GH/GIG routines and an interpolated inverse-CDF table.

## SIS allocation objectives

- tail-probability or conditional-excess optimization;
- total MSE;
- total squared relative error;
- maximum error;
- maximum relative error;
- a selected threshold index.

## Deliberate exclusion

`R2X` only writes a data frame to the Windows clipboard. It is not a numerical
algorithm and is not translated.

R package documentation, list/data-frame presentation, timing wrappers, and
`Runuran` object infrastructure are also excluded. Equivalent numerical values
are exposed through typed Fortran structures.
