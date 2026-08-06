# R to Fortran API map

| R routine | Fortran routine | Notes |
|---|---|---|
| `.beta1Spot` | `beta1_spot` | Elemental pure function |
| `.beta2Spot` | `beta2_spot` | Elemental pure function |
| `.beta1Forward` | `beta1_forward` | Elemental pure function |
| `.beta2Forward` | `beta2_forward` | Elemental pure function |
| `.factorBeta1` | `factor_beta1` | Elemental pure function |
| `.factorBeta2` | `factor_beta2` | Elemental pure function |
| `.NS.estimator` | internal `ns_estimator` | Native QR least squares |
| `.NSS.estimator` | internal `nss_estimator` | Native QR least squares |
| `Nelson.Siegel` | `nelson_siegel_fit` | Rank-1 and rank-2 generic |
| `NSrates` | `ns_rates` | Rank-1 and rank-2 generic |
| `Svensson` | `svensson_fit` | Rank-1 and rank-2 generic |
| `Srates` | `svensson_rates` | `rate_type="spot"` or `"forward"` |
