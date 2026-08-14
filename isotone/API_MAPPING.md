# API mapping

| R `isotone` export | Fortran counterpart | Status |
|---|---|---|
| `gpava` | `gpava_fit`, `gpava_fit_repeated` | Implemented |
| `weighted.median` | `weighted_median_value` | Implemented |
| `weighted.fractile` | `weighted_fractile_value` | Implemented |
| `activeSet` | `active_set`, `active_set_custom` | Implemented |
| `aSolver` | `a_solver`, `ISO_ASYLS` | Implemented |
| `dSolver` | `d_solver`, `ISO_L1` | Implemented |
| `eSolver` | `e_solver`, `ISO_L1EPS` | Implemented |
| `fSolver` | `f_solver`, `active_set_custom` | Implemented with procedure callback/BFGS |
| `hSolver` | `h_solver`, `ISO_HUBER` | Implemented |
| `iSolver` | `i_solver`, `ISO_SILF` | Implemented |
| `lfSolver` | `lf_solver`, `ISO_GLS` | Implemented |
| `lsSolver` | `ls_solver`, `ISO_LS` | Implemented |
| `mSolver` | `m_solver`, `ISO_CHEBYSHEV` | Implemented |
| `oSolver` | `o_solver`, `ISO_LP` | Implemented |
| `pSolver` | `p_solver`, `ISO_QUANTILE` | Implemented |
| `sSolver` | `s_solver`, `ISO_POISSON` | Implemented |
| `mregnn` | `mregnn` | Implemented |
| `mregnnM` | `mregnn_monotone` | Implemented |
| `mregnnP` | `mregnn_positive` | Implemented |
| `print.gpava`, `plot.gpava` | none | R display/plotting omitted |
| `print.activeset`, `summary.activeset` | result fields/KKT values | R display omitted |

The Fortran active-set interface uses integer solver constants instead of R
character dispatch. The mathematical order convention is unchanged:
`[i,j]` means `x(j) >= x(i)`.
