module grpreg
   use grpreg_kinds, only : dp
   use grpreg_types, only : grpreg_fit_type, grpreg_cv_type, grpreg_selection_type, mfdr_result_type, &
      spline_expansion_type
   use grpreg_api, only : grpreg_fit, gbridge_fit, grpsurv_fit, coef_grpreg, predict_grpreg, &
      predict_grpsurv_link, breslow_baseline, predict_grpsurv_survival, predict_grpsurv_hazard, &
      predict_grpsurv_median, loglik_grpreg, residuals_grpreg, select_grpreg, count_nonzero, &
      count_nonzero_groups, group_norms
   use grpreg_cv, only : cv_grpreg, cv_grpsurv, coef_cv_grpreg, predict_cv_grpreg, auc_cv_grpsurv, &
      summarize_cv_grpreg
   use grpreg_mfdr, only : mfdr_grpreg
   use grpreg_spline, only : expand_spline, predict_spline
   use grpreg_data, only : gen_nonlinear_data
   implicit none
   private
   public :: dp
   public :: grpreg_fit_type, grpreg_cv_type, grpreg_selection_type, mfdr_result_type, spline_expansion_type
   public :: grpreg_fit, gbridge_fit, grpsurv_fit
   public :: coef_grpreg, predict_grpreg, predict_grpsurv_link, breslow_baseline, predict_grpsurv_survival
   public :: predict_grpsurv_hazard, predict_grpsurv_median
   public :: loglik_grpreg, residuals_grpreg, select_grpreg, count_nonzero, count_nonzero_groups, group_norms
   public :: cv_grpreg, cv_grpsurv, coef_cv_grpreg, predict_cv_grpreg, auc_cv_grpsurv, summarize_cv_grpreg
   public :: mfdr_grpreg, expand_spline, predict_spline, gen_nonlinear_data
end module grpreg
