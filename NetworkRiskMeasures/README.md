# NetworkRiskMeasures

Modern free-form Fortran translation of the numerical code in CRAN
NetworkRiskMeasures 0.1.7, by Carlos Cinelli and Thiago Cristiano Silva.
GPL-3.0-only; see LICENSE and NOTICE.md.

## Build and test

Place this directory at the top level of
https://github.com/Beliavsky/Fortran-from-R-packages alongside unchanged
`rfortran-core`, `expm`, and `rfortran-linalg` translations, then run:

```sh
fpm build
fpm test
fpm run --example --all
```

The manifest reuses sibling `rfortran-core` and `expm` through path dependencies.
The latter reuses `rfortran-linalg` and its pinned upstream fortran-lapack dependency.
There is no separately installed BLAS/LAPACK requirement, including with FPM
and gfortran on Windows; there are no system numerical-library link flags.
Only this package is distributed in its ZIP. The Windows toolchain is a
supported design target, not a claim of a Windows-host test run.

Import `network_risk_measures` for the complete public API, shared `dp`, result
types and diagnostic constants. Numerical inputs use `real(dp)`. Dense arrays
and explicit status outputs replace R objects, exceptions and messages.
The example demonstrates the real calling convention; API.md documents all
arguments, shapes, defaults and material differences.

## Translation coverage

Package status: **substantial**. **12 of 12 (100%)** exported computational
R functions are mapped. The fraction measures mapped computational functions,
not complete R compatibility or equivalent random streams.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| risk_matrix | network_risk_base | risk_matrix | substantial |
| impact_matrix | network_risk_base | impact_matrix | substantial |
| vulnerability_matrix | network_risk_base | vulnerability_matrix | substantial |
| criticality | network_risk_base | criticality | substantial |
| communicability_matrix | network_communicability | communicability_matrix | substantial |
| impact_susceptibility | network_communicability | impact_susceptibility | substantial |
| impact_fluidity | network_communicability | impact_fluidity | substantial |
| impact_diffusion | network_communicability | impact_diffusion | substantial |
| contagion | network_contagion | contagion | substantial |
| matrix_estimation | network_reconstruction | matrix_estimation | substantial |
| max_ent | network_reconstruction | max_ent | substantial |
| min_dens | network_reconstruction | min_dens | substantial |

The denominator is the 12 distinct `export(...)` names in the preserved
NAMESPACE. Default/igraph adapters and summary/printing/plot methods are only
registered with `S3method(...)`, not separately exported symbols. Dense default
methods implement their generic's numerical mapping. Numeric summary.contagion
calculations are additionally available via `summarize_contagion` and the
contagion result; they are not counted a second time. No exported computational
function is listed as untranslated.

## Directed network conventions

An assets matrix has A(i,j) equal to vertex i's asset claim on j. Vulnerability
is A(i,j)/buffer(i); impact is its transpose, so impact(i,j) describes loss
transmitted from defaulting i to j. Liabilities inputs are transposed before
normalization. Propagation multiplies the **row** stress increment by impact.
Already-normalized impact/vulnerability inputs need no buffer except optional
buffer validation and monetary contagion summaries. Binary risk matrices
threshold at one; criticality and communicability-derived risk measures cap
impact entries above at one. Contagion does **not** cap the impact matrix.

## Important numerical behavior

- Nonnegative finite financial quantities are required. Negative exposures,
  negative/near-zero buffers, nonfinite values and incompatible shapes return
  status codes; this intentionally tightens R's permissive invalid domains.
- Disconnected and zero-total networks are handled explicitly. An isolated
  vertex has zero unweighted susceptibility and diffusion. Self-exposures are
  allowed for risk/contagion inputs; reconstruction enforces a zero diagonal.
- DebtRank supports repeated propagation and original single-hit behavior.
  Threshold cascades support trigger alpha in (0,1]. The source's **raw final
  stress may exceed one**; only propagation and weighted stress summaries cap
  stress. Monetary losses remain uncapped. Omitted-buffer/weight defaults are
  initialized safely, fixing the source use of `v` before construction.
  Scenario histories and step counts
  preserve that source behavior. Nonconvergence is explicit, not silent.
- Full communicability reuses sibling expm. Finite-degree Taylor terms use a
  factorial-free recurrence. Degree zero is identity, correcting upstream's
  accidental R `1:0` loop. R sparse storage/igraph adapters are omitted.
- Maximum entropy is alternating row/column scaling of an off-diagonal prior.
  Infeasible margins are rejected; boundary-feasible margins can converge
  slowly and return the last iterate with `nr_iterations`.
- Minimum density preserves upstream's linear deficit objective weighted by
  **squared** alpha and delta; it does not substitute a squared-error objective.
  It is a stochastic heuristic, not a global minimum-density guarantee.
  Scalar alpha/delta replace R recycling. Use a local integer seed or supplied
  uniform stream for reproducibility. No R RNG sequence parity is claimed.
  Its singleton link-removal indexing bug is corrected. API.md gives details.

## Validation

`fpm test` covers all 12 mappings against executable R 0.1.7 oracle values,
including asymmetric graph orientation, weighted/full/truncated communicability,
three contagion modes, summaries, reconstruction and controlled-draw minimum
density. Independent analytic directed-chain and zero-graph tests, shape/domain
rejections, disconnected vertices, RNG exhaustion/repeatability, rejection/
acceptance, iteration limits, and deterministic dense graph properties supplement
those references. `reference/generate_reference.R` regenerates the oracle module
using the actual installed upstream package; the final tests need no R runtime.
Its controlled-draw minimum-density oracle changes random-draw entry points
only; a separate ordinary R seed-192 result is also preserved.

Default, checked-debug and release FPM build/test/all-example runs are required
before packaging. Never compile this package with `-ffast-math`, `-Ofast`,
`-ffinite-math-only` or equivalent assumptions: validation uses IEEE semantics.
All maintained routines are pure where allowed by dependencies; the four
communicability-derived entry points cannot be pure because the unchanged
shared matrix-exponential API is impure.
