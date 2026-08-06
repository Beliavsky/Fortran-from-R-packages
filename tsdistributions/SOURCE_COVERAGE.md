# Source coverage

All computational exports in the upstream NAMESPACE are represented:

- `authorized_domain`
- `ddist`, `pdist`, `qdist`, `rdist`
- `d/p/q/rstd`, `d/p/q/rsnorm`, `d/p/q/rsstd`
- `d/p/q/rged`, `d/p/q/rsged`
- `d/p/q/rnig`, `d/p/q/rgh`, `d/p/q/rghst`
- `d/p/q/rjsu`
- raw `d/p/q/rghyp`
- `dskewness`, `dkurtosis`
- `nigtransform`, `ghyptransform`
- `distribution_bounds`, `distribution_modelspec`
- `spd_modelspec`, `dspd`, `pspd`, `qspd`, `rspd`

The computational methods behind estimation, covariance, moments, and profiling are exposed as `estimate_distribution`, `information_covariance`, `distribution_moments`, and `tsprofile`.

R-only print, summary, plotting, formula, S3, data-table, TMB loader, and future/parallel orchestration code is not translated.
