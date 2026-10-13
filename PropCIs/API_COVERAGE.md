# API coverage

Status: substantial. 21 of 21 exported computational R functions (100%) are mapped;
this is not a claim of full R interface or optimizer/RNG compatibility.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `acceptbin` | `propcis_frequentist` | `acceptbin` | substantial |
| `add4ci` | `propcis_frequentist` | `add4ci` | substantial |
| `addz2ci` | `propcis_frequentist` | `addz2ci` | substantial |
| `blakerci` | `propcis_frequentist` | `blakerci` | substantial |
| `diffci.bayes` | `propcis_bayesian` | `diffci_bayes` | substantial |
| `diffci.bayes.hpd` | `propcis_bayesian` | `diffci_bayes_hpd` | substantial |
| `diffpropci.Wald.mp` | `propcis_frequentist` | `diffpropci_wald_mp` | substantial |
| `diffpropci.mp` | `propcis_frequentist` | `diffpropci_mp` | substantial |
| `diffscoreci` | `propcis_frequentist` | `diffscoreci` | substantial |
| `exactci` | `propcis_frequentist` | `exactci` | substantial |
| `limit` | `propcis_frequentist` | `limit` | substantial |
| `midPci` | `propcis_frequentist` | `midpci` | substantial |
| `oddsratioci.mp` | `propcis_frequentist` | `oddsratioci_mp` | substantial |
| `orci.bayes` | `propcis_bayesian` | `orci_bayes` | substantial |
| `orscoreci` | `propcis_frequentist` | `orscoreci` | substantial |
| `riskscoreci` | `propcis_frequentist` | `riskscoreci` | substantial |
| `rrci.bayes` | `propcis_bayesian` | `rrci_bayes` | substantial |
| `scoreci` | `propcis_frequentist` | `scoreci` | substantial |
| `scoreci.mp` | `propcis_frequentist` | `scoreci_mp` | substantial |
| `wald2ci` | `propcis_frequentist` | `wald2ci` | substantial |
| `z2stat` | `propcis_frequentist` | `z2stat` | substantial |

## Per-function notes

- `acceptbin` → `acceptbin`: Binomial acceptance probability; scalar integer counts and explicit validation.
- `add4ci` → `add4ci`: Agresti-Coull add-four interval with point estimate; R htest construction omitted.
- `addz2ci` → `addz2ci`: Confidence-dependent pseudo-count interval with point estimate; R htest construction omitted.
- `blakerci` → `blakerci`: Original tolerance-stepped acceptance inversion; tolerance range 1e-9 to 1 and extreme-tail status documented.
- `diffci.bayes` → `diffci_bayes`: Same beta-posterior equal-tail target using deterministic infinite quadrature and bracketed inversion; simulation starts, nsim, R RNG, and optim trajectories omitted.
- `diffci.bayes.hpd` → `diffci_bayes_hpd`: Shortest connected credible interval using quantile-span optimization; fixes asymmetric-prior reflection and upper clipping defects; not a disconnected HPD set or exact optim trajectory.
- `diffpropci.Wald.mp` → `diffpropci_wald_mp`: Unadjusted paired Wald difference interval and estimate; R htest construction omitted.
- `diffpropci.mp` → `diffpropci_mp`: Adjusted paired difference interval and estimate, clipped to [-1,1]; R htest construction omitted.
- `diffscoreci` → `diffscoreci`: Miettinen-Nurminen corrected score inversion; stable constrained likelihood root replaces the cubic helper, preserving original endpoint-search tolerances.
- `exactci` → `exactci`: Clopper-Pearson interval using equivalent beta quantiles; reports limited numerical resolution for extreme confidence tails.
- `limit` → `limit`: Exported multiplicative odds-ratio grid helper; stable quadratic root handles initial ratio one, with bounded iterations and NaN failure.
- `midPci` → `midpci`: Original 0.0001 + k*0.0005 grid and strict comparisons, including analytic zero/all-success branches; extreme-tail resolution flagged.
- `oddsratioci.mp` → `oddsratioci_mp`: Paired odds-ratio score interval; explicitly corrects zero-discordance endpoint roundoff.
- `orci.bayes` → `orci_bayes`: Deterministic beta-logit convolution and equal-tail inversion, equivalent to the F-ratio posterior; fixes upstream all-success undefined-variable and asymmetric-prior complement defects; simulation starts omitted.
- `orscoreci` → `orscoreci`: Original 1.001 multiplicative score grids with stable constrained root and exact zero/infinite support endpoints.
- `riskscoreci` → `riskscoreci`: Same constrained-binomial risk-ratio score using bracketed inversion; replaces unstable cubic and boundary grid searches, so boundary bounds can differ by their former 1.0001 step.
- `rrci.bayes` → `rrci_bayes`: Beta-posterior risk-ratio equal-tail interval from deterministic quadrature and log-ratio inversion; simulation starts, nsim, RNG, and optim trajectories omitted.
- `scoreci` → `scoreci`: Wilson score interval rounded to four decimals; uses Fortran ANINT and a tail-stable normal critical value.
- `scoreci.mp` → `scoreci_mp`: Paired score inversion with original stepping and stopping tolerances; degenerate variance handled explicitly.
- `wald2ci` → `wald2ci`: Independent Wald or add-two-per-sample interval and estimate; fixes upstream unusable vector default by choosing AC.
- `z2stat` → `z2stat`: Exported Miettinen-Nurminen statistic; maximizes the constrained binomial likelihood by a bracketed derivative root instead of the fragile cubic.
