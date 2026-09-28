module car
   use car_anova, only: anova_lm_type3
   use car_kinds, only: dp
   use car_box_tidwell, only: box_tidwell
   use car_contrasts, only: helmert_contrast, sum_contrast, treatment_contrast
   use car_data, only: compare_coefs_numeric, match_coefficient_names, recode_numeric, some_indices, &
                       strings2factors_codes, tapply_numeric, which_names
   use car_diagnostics, only: durbin_watson, levene_test, ncv_test, outlier_bonferroni
   use car_geometry, only: ellipse_coordinates
   use car_hypothesis, only: make_hypothesis
   use car_kernels, only: adaptive_kernel, biweight_kernel, bw_nrd0, epanechnikov_kernel, normal_kernel
   use car_lm, only: confint_lm, lm_bootstrap_case, lm_refit_rows, lm_summary, predict_lm
   use car_manova, only: manova_test_statistics, manova_type3_ssp
   use car_matrix, only: delta_method_gradient, gvif_from_covariance, hccm_matrix, linear_hypothesis_test, &
                         sigma_hat, wcrossprod, wcrossprod_matrix_weight, wcrossprod_vector_weight
   use car_po_test, only: po_test
   use car_power_fit, only: box_cox_profile, inv_transform_estimate, power_transform_fit, power_transform_test
   use car_transforms, only: basic_power, bc_power, bcn_power, bcn_power_inverse, box_cox_variable, &
                             logit_transform, yj_power
   implicit none
   private

   public :: adaptive_kernel
   public :: anova_lm_type3
   public :: basic_power
   public :: bc_power
   public :: bcn_power
   public :: bcn_power_inverse
   public :: biweight_kernel
   public :: box_tidwell
   public :: box_cox_profile
   public :: box_cox_variable
   public :: confint_lm
   public :: compare_coefs_numeric
   public :: bw_nrd0
   public :: delta_method_gradient
   public :: dp
   public :: durbin_watson
   public :: ellipse_coordinates
   public :: epanechnikov_kernel
   public :: gvif_from_covariance
   public :: hccm_matrix
   public :: helmert_contrast
   public :: inv_transform_estimate
   public :: levene_test
   public :: linear_hypothesis_test
   public :: lm_bootstrap_case
   public :: lm_refit_rows
   public :: lm_summary
   public :: logit_transform
   public :: make_hypothesis
   public :: match_coefficient_names
   public :: manova_test_statistics
   public :: manova_type3_ssp
   public :: ncv_test
   public :: normal_kernel
   public :: outlier_bonferroni
   public :: po_test
   public :: power_transform_fit
   public :: power_transform_test
   public :: predict_lm
   public :: recode_numeric
   public :: sigma_hat
   public :: some_indices
   public :: strings2factors_codes
   public :: sum_contrast
   public :: tapply_numeric
   public :: treatment_contrast
   public :: wcrossprod
   public :: wcrossprod_matrix_weight
   public :: wcrossprod_vector_weight
   public :: which_names
   public :: yj_power

end module car
