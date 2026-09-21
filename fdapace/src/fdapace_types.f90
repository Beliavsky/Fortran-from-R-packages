module fdapace_types
    use fdapace_kinds, only : dp
    implicit none
    private

    type, public :: real_vector
        real(dp), allocatable :: v(:)
    end type real_vector

    type, public :: fpca_inputs
        integer, allocatable :: lid(:)
        type(real_vector), allocatable :: ly(:)
        type(real_vector), allocatable :: lt(:)
    end type fpca_inputs

    type, public :: sparse_sample
        type(real_vector), allocatable :: ly(:)
        type(real_vector), allocatable :: lt(:)
    end type sparse_sample

    type, public :: gp_functional_data
        real(dp), allocatable :: y(:,:)
        real(dp), allocatable :: yn(:,:)
        real(dp), allocatable :: phi(:,:)
        real(dp), allocatable :: xi(:,:)
        real(dp), allocatable :: pts(:)
        logical :: has_noisy = .false.
    end type gp_functional_data

    type, public :: sparse_gp_result
        type(sparse_sample) :: sample
        type(sparse_sample) :: true_sample
        real(dp), allocatable :: xi(:,:)
        integer, allocatable :: ni(:)
        logical :: has_true = .false.
    end type sparse_gp_result

    type, public :: fpca_result
        real(dp), allocatable :: lambda(:)
        real(dp), allocatable :: phi(:,:)
        real(dp), allocatable :: xi_est(:,:)
        real(dp), allocatable :: obs_grid(:)
        real(dp), allocatable :: work_grid(:)
        real(dp), allocatable :: mu(:)
        real(dp), allocatable :: smoothed_cov(:,:)
        real(dp), allocatable :: fitted_cov(:,:)
        real(dp), allocatable :: fitted_corr(:,:)
        real(dp), allocatable :: cum_fve(:)
        real(dp), allocatable :: fitted_y(:,:)
        real(dp) :: sigma2 = 0.0_dp
        real(dp) :: fve = 0.0_dp
        integer :: select_k = 0
        logical :: has_sigma2 = .false.
    end type fpca_result

    type, public :: fsvd_result
        real(dp), allocatable :: cr_cov(:,:)
        real(dp), allocatable :: s_values(:)
        real(dp), allocatable :: can_corr(:)
        real(dp), allocatable :: s_fun1(:,:)
        real(dp), allocatable :: grid1(:)
        real(dp), allocatable :: scores1(:,:)
        real(dp), allocatable :: s_fun2(:,:)
        real(dp), allocatable :: grid2(:)
        real(dp), allocatable :: scores2(:,:)
        real(dp) :: fve = 0.0_dp
        integer :: nsvd = 0
    end type fsvd_result

    type, public :: select_k_result
        integer :: k = 0
        real(dp) :: criterion = 0.0_dp
        logical :: has_criterion = .false.
    end type select_k_result

    type, public :: bwnn_result
        real(dp) :: cov_bw = 0.0_dp
        real(dp) :: mu_bw = 0.0_dp
        logical :: has_cov = .false.
        logical :: has_mean = .false.
    end type bwnn_result

    type, public :: fccor_result
        real(dp), allocatable :: corr(:)
        real(dp), allocatable :: tout(:)
        real(dp) :: bw(5) = 0.0_dp
    end type fccor_result

    type, public :: dyn_test_result
        real(dp) :: stats = 0.0_dp
        real(dp) :: pval = 0.0_dp
    end type dyn_test_result

    type, public :: mean_curve_result
        real(dp), allocatable :: mu(:)
        real(dp), allocatable :: work_grid(:)
        real(dp) :: bw_mu = 0.0_dp
    end type mean_curve_result

    type, public :: cov_surface_result
        real(dp), allocatable :: cov(:,:)
        real(dp), allocatable :: work_grid(:)
        real(dp) :: sigma2 = 0.0_dp
        real(dp) :: bw_cov = 0.0_dp
    end type cov_surface_result

    type, public :: mean_ci_result
        real(dp), allocatable :: grid(:)
        real(dp), allocatable :: lower(:)
        real(dp), allocatable :: upper(:)
        real(dp) :: level = 0.95_dp
    end type mean_ci_result

    type, public :: fpca_der_result
        real(dp), allocatable :: mu_der(:)
        real(dp), allocatable :: phi_der(:,:)
        integer :: derivative_order = 1
    end type fpca_der_result

    type, public :: fvpa_result
        type(fpca_result) :: fpca_y
        type(fpca_result) :: fpca_r
        real(dp) :: sigma2 = 0.0_dp
        real(dp) :: delta = 0.0_dp
    end type fvpa_result

    type, public :: stringing_result
        integer, allocatable :: order(:)
        real(dp), allocatable :: stringed_x(:,:)
        real(dp), allocatable :: standardized_x(:,:)
        real(dp), allocatable :: distance(:,:)
        logical :: standardized = .false.
    end type stringing_result

    type, public :: cluster_result
        integer, allocatable :: cluster(:)
        real(dp), allocatable :: centers(:,:)
        integer :: iterations = 0
        logical :: converged = .false.
    end type cluster_result

    type, public :: fopt_result
        integer, allocatable :: indices(:)
        real(dp), allocatable :: opt_des(:)
        real(dp) :: r2 = 0.0_dp
        real(dp) :: r2_adj = 0.0_dp
        real(dp) :: ridge = 0.0_dp
    end type fopt_result

    type, public :: sbf_result
        real(dp), allocatable :: fit(:,:)
        real(dp), allocatable :: nw(:,:)
        real(dp) :: mean_y = 0.0_dp
        integer :: iterations = 0
        real(dp) :: error = 0.0_dp
        logical :: converged = .false.
    end type sbf_result

    type, public :: fam_result
        real(dp), allocatable :: fam(:,:)
        real(dp), allocatable :: xi(:,:)
        real(dp), allocatable :: bw(:)
        real(dp), allocatable :: lambda(:)
        real(dp), allocatable :: phi(:,:)
        real(dp), allocatable :: work_grid(:)
        real(dp) :: mu = 0.0_dp
    end type fam_result

    type, public :: multifam_result
        real(dp), allocatable :: sbfit(:,:)
        real(dp), allocatable :: xi(:,:)
        real(dp), allocatable :: bw(:)
        real(dp), allocatable :: lambda(:)
        real(dp) :: mu = 0.0_dp
    end type multifam_result

    type, public :: flm_result
        real(dp), allocatable :: alpha(:)
        real(dp), allocatable :: beta(:,:)
        real(dp), allocatable :: y_hat(:,:)
        real(dp), allocatable :: y_pred(:,:)
        real(dp) :: r2 = 0.0_dp
        real(dp) :: p_value = 0.0_dp
        logical :: has_p_value = .false.
    end type flm_result

    type, public :: flm_ci_result
        real(dp), allocatable :: alpha_lower(:)
        real(dp), allocatable :: alpha_upper(:)
        real(dp), allocatable :: beta_lower(:,:)
        real(dp), allocatable :: beta_upper(:,:)
        real(dp) :: level = 0.95_dp
    end type flm_ci_result

    type, public :: fpc_quantile_result
        real(dp), allocatable :: pred_quantile(:,:)
        real(dp), allocatable :: pred_cdf(:,:)
        real(dp), allocatable :: beta(:,:)
        real(dp), allocatable :: cdf_grid(:)
        real(dp), allocatable :: quantiles(:)
    end type fpc_quantile_result

    type, public :: fcreg_result
        real(dp), allocatable :: beta(:,:)
        real(dp), allocatable :: beta0(:)
        real(dp), allocatable :: r2(:)
        real(dp), allocatable :: out_grid(:)
    end type fcreg_result

    type, public :: tvam_result
        real(dp), allocatable :: mean_t(:)
        real(dp), allocatable :: components(:,:,:)
        real(dp), allocatable :: grid_t(:)
        real(dp), allocatable :: x_eval(:,:)
    end type tvam_result

    type, public :: vcam_result
        real(dp), allocatable :: fitted_y(:,:)
        real(dp), allocatable :: phi_est(:,:)
        real(dp), allocatable :: beta0_est(:)
        real(dp), allocatable :: beta_est(:,:)
        real(dp), allocatable :: grid_t(:)
        real(dp), allocatable :: grid_x(:,:)
    end type vcam_result

    type, public :: wfda_result
        real(dp), allocatable :: h(:,:)
        real(dp), allocatable :: h_inv(:,:)
        real(dp), allocatable :: aligned(:,:)
        real(dp), allocatable :: costs(:)
        real(dp) :: lambda = 0.0_dp
    end type wfda_result
end module fdapace_types
