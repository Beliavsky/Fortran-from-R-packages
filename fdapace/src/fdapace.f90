module fdapace
    use fdapace_kinds, only : dp
    use fdapace_types, only : bwnn_result, cluster_result, cov_surface_result, dyn_test_result, fam_result, &
        fccor_result, fcreg_result, flm_ci_result, flm_result, fopt_result, fpc_quantile_result, fpca_der_result, &
        fpca_inputs, fpca_result, fsvd_result, fvpa_result, gp_functional_data, mean_ci_result, mean_curve_result, &
        multifam_result, real_vector, sbf_result, select_k_result, sparse_gp_result, sparse_sample, stringing_result, &
        tvam_result, vcam_result, wfda_result
    use fdapace_basis, only : create_basis
    use fdapace_smoothing, only : convert_support, convert_support_covariance, convert_support_functions, &
        cumtrapz_rcpp, lwls1d, lwls2d, lwls2d_deriv, norm_curv_to_area, trapz_rcpp
    use fdapace_statistics, only : bw_nn, dyn_corr, dyn_test, fc_cor
    use fdapace_growth, only : make_bw_to_zscore_02y, make_hc_to_zscore_02y, make_ln_to_zscore_02y
    use fdapace_simulation, only : make_fpca_inputs_dense, make_gp_functional_data, make_sparse_gp, sparsify, wiener
    use fdapace_fpca, only : fpca_dense, fsvd_dense, get_normalised_sample, &
        get_normalized_sample, select_k_fixed, select_k_fve
    use fdapace_covariance, only : get_cov_surface, get_cr_cor_yx, get_cr_cor_yz, get_cr_cov_yx, get_cr_cov_yz, &
        get_mean_ci, get_mean_curve
    use fdapace_analysis, only : fclust, fpca_der, fvpa, k_cfc, stringing, wfda
    use fdapace_models, only : fam, fc_reg, flm_ci_scalar_dense, flm_scalar_dense, f_opt_des, fpc_quantile, &
        multi_fam, sb_fitting, tvam, vcam
    use fdapace_math, only : seed_rng
    implicit none
    private

    public :: bwnn_result
    public :: bw_nn
    public :: cluster_result
    public :: convert_support
    public :: convert_support_covariance
    public :: convert_support_functions
    public :: cov_surface_result
    public :: create_basis
    public :: cumtrapz_rcpp
    public :: dp
    public :: dyn_corr
    public :: dyn_test
    public :: dyn_test_result
    public :: fam
    public :: fam_result
    public :: fc_cor
    public :: fccor_result
    public :: fc_reg
    public :: fclust
    public :: fcreg_result
    public :: flm_ci_result
    public :: flm_ci_scalar_dense
    public :: flm_result
    public :: flm_scalar_dense
    public :: f_opt_des
    public :: fopt_result
    public :: fpc_quantile
    public :: fpc_quantile_result
    public :: fpca_dense
    public :: fpca_der
    public :: fpca_der_result
    public :: fpca_inputs
    public :: fpca_result
    public :: fsvd_dense
    public :: fsvd_result
    public :: fvpa
    public :: fvpa_result
    public :: get_cov_surface
    public :: get_cr_cor_yx
    public :: get_cr_cor_yz
    public :: get_cr_cov_yx
    public :: get_cr_cov_yz
    public :: get_mean_ci
    public :: get_mean_curve
    public :: get_normalised_sample
    public :: get_normalized_sample
    public :: gp_functional_data
    public :: k_cfc
    public :: lwls1d
    public :: lwls2d
    public :: lwls2d_deriv
    public :: make_bw_to_zscore_02y
    public :: make_fpca_inputs_dense
    public :: make_gp_functional_data
    public :: make_hc_to_zscore_02y
    public :: make_ln_to_zscore_02y
    public :: make_sparse_gp
    public :: mean_ci_result
    public :: mean_curve_result
    public :: multi_fam
    public :: multifam_result
    public :: norm_curv_to_area
    public :: real_vector
    public :: sb_fitting
    public :: sbf_result
    public :: seed_rng
    public :: select_k_fixed
    public :: select_k_fve
    public :: select_k_result
    public :: sparse_gp_result
    public :: sparse_sample
    public :: sparsify
    public :: stringing
    public :: stringing_result
    public :: trapz_rcpp
    public :: tvam
    public :: tvam_result
    public :: vcam
    public :: vcam_result
    public :: wfda
    public :: wfda_result
    public :: wiener
end module fdapace
