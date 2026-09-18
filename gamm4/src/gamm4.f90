module gamm4
   use gamm4_kinds, only : dp
   use gamm4_types, only : gamm4_family_gaussian, gamm4_smooth_t, gamm4_control_t, &
      gamm4_result_t, gamm4_vb_result_t
   use gamm4_covariance, only : gamm4_get_vb
   use gamm4_fit_mod, only : gamm4_fit
   use lme4, only : random_term_t, family_binomial, family_poisson, family_gamma, &
      family_inverse_gaussian, family_negative_binomial, covariance_unstructured, &
      covariance_diagonal, covariance_compound_symmetry, covariance_ar1
   implicit none
   public
end module gamm4
