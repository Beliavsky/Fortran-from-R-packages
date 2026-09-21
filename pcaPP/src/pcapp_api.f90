module pcapp_api
  use pcapp_kinds, only: dp
  use pcapp_types, only: scale_result, median_result, pca_result, covariance_result, tuning_result
  use pcapp_stats, only: qn_scale, cor_fk_vector, cor_fk_matrix
  use pcapp_l1median, only: l1median, l1median_bfgs, l1median_cg, l1median_hocr
  use pcapp_l1median, only: l1median_nlm, l1median_nm, l1median_vazh
  use pcapp_pca, only: scale_adv, pca_grid, spca_grid, pca_proj
  use pcapp_pca, only: cov_pc, cov_pca_grid, cov_pca_proj, opt_tpo, opt_bic, data_zou
  implicit none
  private

  public :: dp
  public :: scale_result, median_result, pca_result, covariance_result, tuning_result
  public :: qn, cor_fk, cor_fk_vec_api, cor_fk_mat_api
  public :: scale_adv, pca_grid, spca_grid, pca_proj
  public :: cov_pc, cov_pca_grid, cov_pca_proj, data_zou
  public :: l1median, l1median_bfgs, l1median_cg, l1median_hocr
  public :: l1median_nlm, l1median_nm, l1median_vazh
  public :: opt_tpo, opt_bic

  interface cor_fk
    module procedure cor_fk_vec_api
    module procedure cor_fk_mat_api
  end interface cor_fk

contains

  pure real(dp) function qn(x, corr_fact) result(value)
    real(dp), intent(in) :: x(:) !! Numeric sample for the pcaPP Qn robust scale estimator.
    real(dp), intent(in), optional :: corr_fact !! Optional asymptotic normalization constant replacing the default.

    if (present(corr_fact)) then
      value = qn_scale(x, corr_fact)
    else
      value = qn_scale(x)
    end if
  end function qn

  pure real(dp) function cor_fk_vec_api(x, y) result(value)
    real(dp), intent(in) :: x(:) !! First vector for the pcaPP Kendall tau-b correlation.
    real(dp), intent(in) :: y(:) !! Second vector, with the same length as x.

    value = cor_fk_vector(x, y)
  end function cor_fk_vec_api

  pure function cor_fk_mat_api(x) result(corr)
    real(dp), intent(in) :: x(:,:) !! Observation-by-variable matrix for pairwise pcaPP Kendall correlations.
    real(dp) :: corr(size(x, 2), size(x, 2))

    call cor_fk_matrix(x, corr)
  end function cor_fk_mat_api

end module pcapp_api
