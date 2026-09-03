! Public facade for the modern Fortran translation of geepack 1.3-13.
! Upstream license: GPL (>= 3). See LICENSE, NOTICE.md, and PROVENANCE.md.
! Exact upstream copyright-holder attribution is preserved in NOTICE.md and upstream/DESCRIPTION.
module geepack
   use r_kinds, only : dp
   use geepack_status, only : GEE_OK, GEE_ERR_SHAPE, GEE_ERR_ARGUMENT, GEE_ERR_SINGULAR, &
      GEE_ERR_MAXITER, GEE_ERR_INVALID_MEAN, GEE_ERR_CORRELATION
   use geepack_links, only : LINK_IDENTITY, LINK_LOGIT, LINK_PROBIT, LINK_CLOGLOG, LINK_LOG, &
      LINK_RECIPROCAL, LINK_FISHERZ, LINK_LWYBC2, LINK_LWYLOG, VAR_GAUSSIAN, VAR_BINOMIAL, &
      VAR_POISSON, VAR_GAMMA, link_function, link_inverse, link_derivative, variance_function, &
      variance_derivative, valid_mean
   use geepack_correlations, only : COR_INDEPENDENCE, COR_EXCHANGEABLE, COR_AR1, COR_UNSTRUCTURED, &
      COR_USERDEFINED, COR_FIXED, working_correlation, correlation_rho_derivative
   use geepack_design, only : gen_zcor, gen_zodds, fixed_to_zcor
   use geepack_gee, only : gee_spec, gee_result, fit_geese
   use geepack_ordinal, only : ordinal_spec, ordinal_result, fit_ordgee, odds_to_p11, &
      p11_odds_derivative, p11_mean_derivatives
   use geepack_metrics, only : qic_result, quasi_likelihood, compute_qic, compare_coefficients, fitted_means
   use geepack_relative_risk, only : make_relative_risk_copy, fit_relative_risk
   use geepack_inference, only : coefficient_wald_summary, wald_contrast, chi_square_survival
   implicit none
   private

   public :: dp
   public :: GEE_OK, GEE_ERR_SHAPE, GEE_ERR_ARGUMENT, GEE_ERR_SINGULAR
   public :: GEE_ERR_MAXITER, GEE_ERR_INVALID_MEAN, GEE_ERR_CORRELATION
   public :: LINK_IDENTITY, LINK_LOGIT, LINK_PROBIT, LINK_CLOGLOG, LINK_LOG
   public :: LINK_RECIPROCAL, LINK_FISHERZ, LINK_LWYBC2, LINK_LWYLOG
   public :: VAR_GAUSSIAN, VAR_BINOMIAL, VAR_POISSON, VAR_GAMMA
   public :: COR_INDEPENDENCE, COR_EXCHANGEABLE, COR_AR1, COR_UNSTRUCTURED
   public :: COR_USERDEFINED, COR_FIXED
   public :: link_function, link_inverse, link_derivative, variance_function, variance_derivative, valid_mean
   public :: working_correlation, correlation_rho_derivative
   public :: gen_zcor, gen_zodds, fixed_to_zcor
   public :: gee_spec, gee_result, fit_geese
   public :: ordinal_spec, ordinal_result, fit_ordgee
   public :: odds_to_p11, p11_odds_derivative, p11_mean_derivatives
   public :: qic_result, quasi_likelihood, compute_qic, compare_coefficients, fitted_means
   public :: make_relative_risk_copy, fit_relative_risk
   public :: coefficient_wald_summary, wald_contrast, chi_square_survival

end module geepack
