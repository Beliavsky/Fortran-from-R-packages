# GPareto R -> Fortran API mapping

| Upstream GPareto routine | Fortran counterpart | Status |
|---|---|---|
| `CPF` | `compute_cpf`, `cpf_result` | Implemented (computation only) |
| `ZDT1`, `ZDT2`, `ZDT3`, `ZDT4`, `ZDT6` | `zdt1` ... `zdt6` | Implemented |
| `P1`, `P2` | `p1_test`, `p2_test` | Implemented |
| `MOP2`, `MOP3` | `mop2`, `mop3` | Implemented |
| `DTLZ1`, `DTLZ2`, `DTLZ3`, `DTLZ7` | `dtlz1` ... `dtlz7` | Implemented |
| `OKA1` | `oka1` | Implemented |
| `nonDomSet` / native `nonDomInd_cpp` | `nondom_set`, `nondominated_indices` | Implemented |
| native `distcpp_2` | `squared_distances` | Implemented |
| `predict_kms` | `predict_gps` | Implemented |
| `checkPredict` | `check_predict` | Euclidean-distance guard implemented |
| `crit_EHI` | `crit_ehi` | Exact 2D; SAA for >2 objectives |
| `crit_qEHI` | `crit_qehi` | Implemented with full batch covariance |
| `crit_EMI` | `crit_emi` | SAA implementation; semi-analytic 2D branch not ported |
| `crit_SMS` | `crit_sms` | Implemented |
| `crit_SUR` | `crit_sur` | Implemented by direct posterior Monte Carlo SUR |
| `crit_optimizer` | `crit_optimizer` | Differential-evolution search |
| `GParetoptim` | `gparetoptim` | Implemented |
| `easyGParetoptim` | `easy_gparetoptim` | Implemented |
| `getDesign` | `get_design` | Implemented |
| `integration_design_optim` | `integration_design_optim` | Halton/MC/SUR importance modes |
| `ParetoSetDensity` | `pareto_set_density` | Native conditional simulation + KDE |
| `prob.of.non.domination` | `probability_nondomination` | Exact 2D/3D |
| `VorobThreshold`, `VorobExpect`, `VorobDev` | `vorob_threshold`, `vorob_expectation`, `vorob_deviation` | Implemented |
| `fastfun` | none in v0.1.0 | Future: mixed deterministic/GP objective adapter |
| plotting functions | none | Intentionally omitted |

## Differences worth noting

- Upstream's special semi-analytical bi-objective EMI formula is not used;
  `crit_emi` uses the same SAA/maximin definition in every dimension.
- Upstream SUR contains specialized closed/semi-closed 2D expected-excursion
  calculations. The Fortran routine evaluates the same one-step uncertainty
  reduction target by direct Gaussian conditioning and Monte Carlo over the
  prospective observation.
- `checkPredict` currently implements the Euclidean guard. The R package's
  covariance-distance option is not separately exposed.
- Upstream `sobol` integration designs map to the native Halton generator in
  this release. This affects point sequence reproducibility, not the target
  integral.
- Optimizer trajectories are not expected to match R because `genoud` and PSO
  are replaced by differential evolution.
