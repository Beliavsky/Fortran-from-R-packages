# Translation coverage

| Upstream export | Fortran routine | Status |
|---|---|---|
| `estimate_sample_moments` | `estimate_sample_moments` | Translated |
| `estimate_skew_t` | `estimate_skew_t` | Translated using adapted fitHeavyTail code |
| `eval_portfolio_moments` | generic `eval_portfolio_moments` | Translated for both parameter types |
| `design_MVSK_portfolio_via_sample_moments` | `design_mvsk_portfolio_via_sample_moments` | All 3 methods represented |
| `design_MVSK_portfolio_via_skew_t` | `design_mvsk_portfolio_via_skew_t` | All 6 methods represented |
| `design_MVSKtilting_portfolio_via_sample_moments` | `design_mvsktilting_portfolio_via_sample_moments` | Both methods represented |

Internal translated/adapted functionality includes simplex projection, PSD
Hessian approximation, convex simplex QP solution, sample moment gradients and
Hessians, skew-t moment gradients and Hessians, SQUAREM/RFPA acceleration,
tracking-error projection, and the skew-t EM/PX-EM numerical dependency.
