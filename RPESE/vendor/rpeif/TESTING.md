# Testing

Five standalone test programs are included.

1. `test_nuisance_stats`
   - normal inverse/CDF consistency
   - nuisance identities
   - type-7 quantiles
   - sample SD
   - lower and upper partial moments
2. `test_influence_core`
   - mean, SD, semi-SD, LPM, Omega, Sharpe, Sortino, and downside Sharpe
   - finite outputs
   - source-compatible versus corrected Omega and Sharpe paths
3. `test_tail_ratios`
   - VaR, ES, ES ratio, VaR ratio, Rachev ratio, and robust mean
   - nuisance-based Rachev shape
   - empirical VaR-ratio compatibility switch
4. `test_robust_prewhiten`
   - robust location and scale
   - outlier clipping
   - AR(1) coefficient and residual calculation
5. `test_dispatch_shape`
   - series dispatcher
   - cleaning and prewhitening integration
   - automatic shape grid
   - estimator aliases and invalid names

Run:

```text
fpm test
```

or:

```text
make MODE=checked test
make MODE=optimized test
```
