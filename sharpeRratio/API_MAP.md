# API map

| Upstream R/C++ routine | Fortran routine or type | Status |
|---|---|---|
| `num_records_up` | `num_records_up` | Direct translation |
| `num_records_down` | `num_records_down` | Direct translation |
| `computeR0bar` | `compute_r0bar`, `computeR0bar` | Direct translation with deterministic seed option |
| `Test.N` | `test_n` | Direct translation using R type-7 quantiles and sample variance |
| `a` | `calibration_a` | Exact serialized spline coefficients |
| `a_medium` | `calibration_a_medium` | Exact serialized spline coefficients |
| `f` | `calibration_f` | Exact serialized spline coefficients |
| `a_full` | `a_full` | Direct translation |
| `f_full` | `f_full` | Direct translation |
| `b` | `correction_b` | Direct translation |
| `theta` | `theta_snr` | Direct translation |
| `ghyp::fit.tuv` | `estimate_tail_exponent` using `fit_ghyp_uv(...,"student")` | Same model family; optimizer supplied by ghyp-fortran |
| `estimateSNR` | `estimate_snr`, `estimateSNR` | Direct numerical translation |
| returned R list | `snr_result` | Typed Fortran replacement |
| `computeR0bar` returned list | `r0_result` | Typed Fortran replacement |
