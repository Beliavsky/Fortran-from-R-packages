# API map

`GenBinomApps` 1.2.1 exports eight computational functions. All are translated.

| R API | Fortran API |
|---|---|
| `dgbinom` | `dgbinom`, `dgbinom_vec` |
| `pgbinom` | `pgbinom`, `pgbinom_vec` |
| `qgbinom` | `qgbinom`, `qgbinom_vec` |
| `rgbinom` | `rgbinom` |
| `clopper.pearson.ci` | `clopper_pearson_ci` |
| `cm.clopper.pearson.ci` | `cm_clopper_pearson_ci` |
| `n.clopper.pearson` | `n_clopper_pearson` |
| `cm.n.clopper.pearson` | `cm_n_clopper_pearson` |

Confidence intervals are returned as the derived type `confidence_interval`.
Sample-size routines return `integer(int64)` values.

The internal beta CDF/quantile and scalar root solver are implemented locally;
there are no runtime package dependencies.
