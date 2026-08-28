# API map

Upstream package: `nleqslv` 3.3.7.

| R API | Fortran API | Status |
|---|---|---|
| `nleqslv` | `solve_nleqslv` | translated |
| `searchZeros` | `search_zeros` | translated |
| `testnslv` | `test_nleqslv` | computational portion translated |
| `print.test.nleqslv` | -- | omitted; presentation only |

## Solver options

R strings are represented by integer constants:

| R value | Fortran constant |
|---|---|
| `method="Newton"` | `NLEQ_NEWTON` |
| `method="Broyden"` | `NLEQ_BROYDEN` |
| `global="none"` | `NLEQ_NONE` |
| `global="cline"` | `NLEQ_CLINE` |
| `global="qline"` | `NLEQ_QLINE` |
| `global="gline"` | `NLEQ_GLINE` |
| `global="dbldog"` | `NLEQ_DBLDOG` |
| `global="pwldog"` | `NLEQ_PWLDOG` |
| `global="hook"` | `NLEQ_HOOK` |
| `xscalm="fixed"` | `NLEQ_SCALE_FIXED` |
| `xscalm="auto"` | `NLEQ_SCALE_AUTO` |

`nleq_options` carries `ftol`, `xtol`, `btol`, `cndtol`, `stepmax`, `delta`,
`sigma`, `maxit`, scaling, band widths, singular-Jacobian handling, Jacobian
checking, tracing, and final-Jacobian return selection.

`delta=-1` corresponds to the upstream `"cauchy"` initialization and
`delta=-2` to `"newton"`.

## Results

`nleq_result` contains the final `x`, `fvec`, termination code/message,
scaling vector, function/Jacobian counts, iteration count, inverse-condition
estimate, and optionally the final Jacobian/Broyden matrix.
