# API coverage

## Translation coverage

Package status: **substantial**. **42 of 42 (100%)** exported computational R functions have numerical mappings. The fraction measures mapped computational functions, not complete R compatibility. Forty mappings are substantial; the two cosinor mappings are partial because rank-deficient coefficient selection differs from R.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `row_t_equalvar` | `matrixtests` | `row_t_equalvar` | substantial |
| `col_t_equalvar` | `matrixtests` | `col_t_equalvar` | substantial |
| `row_t_paired` | `matrixtests` | `row_t_paired` | substantial |
| `col_t_paired` | `matrixtests` | `col_t_paired` | substantial |
| `row_t_onesample` | `matrixtests` | `row_t_onesample` | substantial |
| `col_t_onesample` | `matrixtests` | `col_t_onesample` | substantial |
| `row_t_welch` | `matrixtests` | `row_t_welch` | substantial |
| `col_t_welch` | `matrixtests` | `col_t_welch` | substantial |
| `row_cor_pearson` | `matrixtests` | `row_cor_pearson` | substantial |
| `col_cor_pearson` | `matrixtests` | `col_cor_pearson` | substantial |
| `row_f_var` | `matrixtests` | `row_f_var` | substantial |
| `col_f_var` | `matrixtests` | `col_f_var` | substantial |
| `row_oneway_equalvar` | `matrixtests` | `row_oneway_equalvar` | substantial |
| `col_oneway_equalvar` | `matrixtests` | `col_oneway_equalvar` | substantial |
| `row_oneway_welch` | `matrixtests` | `row_oneway_welch` | substantial |
| `col_oneway_welch` | `matrixtests` | `col_oneway_welch` | substantial |
| `row_kruskalwallis` | `matrixtests` | `row_kruskalwallis` | substantial |
| `col_kruskalwallis` | `matrixtests` | `col_kruskalwallis` | substantial |
| `row_waerden` | `matrixtests` | `row_waerden` | substantial |
| `col_waerden` | `matrixtests` | `col_waerden` | substantial |
| `row_bartlett` | `matrixtests` | `row_bartlett` | substantial |
| `col_bartlett` | `matrixtests` | `col_bartlett` | substantial |
| `row_jarquebera` | `matrixtests` | `row_jarquebera` | substantial |
| `col_jarquebera` | `matrixtests` | `col_jarquebera` | substantial |
| `row_flignerkilleen` | `matrixtests` | `row_flignerkilleen` | substantial |
| `col_flignerkilleen` | `matrixtests` | `col_flignerkilleen` | substantial |
| `row_levene` | `matrixtests` | `row_levene` | substantial |
| `col_levene` | `matrixtests` | `col_levene` | substantial |
| `row_brownforsythe` | `matrixtests` | `row_brownforsythe` | substantial |
| `col_brownforsythe` | `matrixtests` | `col_brownforsythe` | substantial |
| `row_wilcoxon_onesample` | `matrixtests` | `row_wilcoxon_onesample` | substantial |
| `col_wilcoxon_onesample` | `matrixtests` | `col_wilcoxon_onesample` | substantial |
| `row_wilcoxon_paired` | `matrixtests` | `row_wilcoxon_paired` | substantial |
| `col_wilcoxon_paired` | `matrixtests` | `col_wilcoxon_paired` | substantial |
| `row_wilcoxon_twosample` | `matrixtests` | `row_wilcoxon_twosample` | substantial |
| `col_wilcoxon_twosample` | `matrixtests` | `col_wilcoxon_twosample` | substantial |
| `row_cosinor` | `matrixtests` | `row_cosinor` | partial |
| `col_cosinor` | `matrixtests` | `col_cosinor` | partial |
| `row_andersondarling` | `matrixtests` | `row_andersondarling` | substantial |
| `col_andersondarling` | `matrixtests` | `col_andersondarling` | substantial |
| `row_kolmogorovsmirnov_twosample` | `matrixtests` | `row_kolmogorovsmirnov_twosample` | substantial |
| `col_kolmogorovsmirnov_twosample` | `matrixtests` | `col_kolmogorovsmirnov_twosample` | substantial |
