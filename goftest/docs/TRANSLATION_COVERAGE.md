# Translation coverage

| Upstream routine | Fortran counterpart | Status |
|---|---|---|
| `pAD` | `p_ad`, `ad_cdf` | translated |
| `qAD` | `q_ad`, `ad_quantile` | translated |
| `pCvM` | `p_cvm`, `cvm_cdf` | translated |
| `qCvM` | `q_cvm`, `cvm_quantile` | translated |
| `ad.test` | `ad_test`, `ad_test_values` | numerical core translated |
| `cvm.test` | `cvm_test`, `cvm_test_values` | numerical core translated |
| `simpleADtest` | `ad_test_uniform` | translated |
| `simpleCvMtest` | `cvm_test_uniform` | translated |
| `braun` | internal Braun path in `goftest_testing` | translated |
| `recogniseCdf` | `recognise_cdf` | translated |
| `getCdf`, `getfunky` | CDF procedure callback | replaced by native Fortran mechanism |

There is no plotting code in the upstream package. R-specific `htest` object
formatting, namespace lookup, expression deparsing, and argument-name reporting
are interface/presentation machinery and are intentionally omitted.
