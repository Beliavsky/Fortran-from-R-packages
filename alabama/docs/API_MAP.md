# API map

| R function | Fortran API | Notes |
|---|---|---|
| `auglag` | `auglag` | Dispatches equality/inequality/mixed constraints |
| `auglag1` | `auglag1` | Equality-only augmented Lagrangian |
| `auglag2` | `auglag2` | Inequality-only augmented Lagrangian |
| `auglag3` | `auglag3` | Equality + inequality augmented Lagrangian |
| `constrOptim.nl` | `constr_optim_nl` | Dispatches legacy algorithms |
| `adpbar` | `adpbar` | Adaptive barrier for inequalities |
| `augpen` | `augpen` | Augmented penalty for equalities |
| `alabama` (internal R helper) | `alabama_legacy` | Combined barrier + equality penalty |

`alabama_outer_control_t` represents outer-loop controls and
`alabama_inner_control_t` represents inner-optimizer controls.
`alabama_result_t` contains parameters, objective value, convergence code,
constraint values, multipliers, penalty parameter, evaluation counts and KKT
information where applicable.
