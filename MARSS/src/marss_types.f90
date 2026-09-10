! SPDX-License-Identifier: GPL-2.0-only
module marss_types
   use marss_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: marss_parameter_name_len = 64

   type, public :: marss_model_spec
      character(len=16) :: form = "marss"
      integer :: nstate = 0
      character(len=32) :: z = "identity"
      character(len=32) :: a = "scaling"
      character(len=32) :: r = "diagonal and equal"
      character(len=32) :: b = "identity"
      character(len=32) :: u = "unconstrained"
      character(len=32) :: q = "diagonal and unequal"
      character(len=32) :: x0 = "unconstrained"
      character(len=32) :: v0 = "zero"
      character(len=32) :: g = "identity"
      character(len=32) :: h = "identity"
      character(len=32) :: l = "identity"
      character(len=32) :: c = "zero"
      character(len=32) :: d = "zero"
      integer :: tinitx = 0
      logical :: diffuse = .false.
   end type marss_model_spec

   type, public :: marss_dfa_spec
      integer :: ntrends = 1
      character(len=32) :: b = "identity"
      character(len=32) :: q = "identity"
      character(len=32) :: a = "zero"
      character(len=32) :: r = "diagonal and equal"
      character(len=32) :: x0 = "zero"
      character(len=32) :: v0 = "fixed5"
      character(len=32) :: d = "auto"
      integer :: tinitx = 0
      logical :: diffuse = .false.
      logical :: demean = .true.
      logical :: z_score = .true.
   end type marss_dfa_spec

   type, public :: marss_model
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: b(:, :)
      real(dp), allocatable :: u(:)
      real(dp), allocatable :: q(:, :)
      real(dp), allocatable :: z(:, :)
      real(dp), allocatable :: a(:)
      real(dp), allocatable :: r(:, :)
      real(dp), allocatable :: x0(:)
      real(dp), allocatable :: v0(:, :)
      real(dp), allocatable :: b_t(:, :, :)
      real(dp), allocatable :: u_t(:, :)
      real(dp), allocatable :: q_t(:, :, :)
      real(dp), allocatable :: z_t(:, :, :)
      real(dp), allocatable :: a_t(:, :)
      real(dp), allocatable :: r_t(:, :, :)
      real(dp), allocatable :: state_covariates(:, :)
      real(dp), allocatable :: obs_covariates(:, :)
      real(dp), allocatable :: c_coef(:, :)
      real(dp), allocatable :: d_coef(:, :)
      real(dp), allocatable :: c_coef_t(:, :, :)
      real(dp), allocatable :: d_coef_t(:, :, :)
      real(dp), allocatable :: g(:, :)
      real(dp), allocatable :: h(:, :)
      real(dp), allocatable :: l(:, :)
      real(dp), allocatable :: g_t(:, :, :)
      real(dp), allocatable :: h_t(:, :, :)
      real(dp), allocatable :: q_noise(:, :)
      real(dp), allocatable :: r_noise(:, :)
      real(dp), allocatable :: v0_noise(:, :)
      real(dp), allocatable :: q_noise_t(:, :, :)
      real(dp), allocatable :: r_noise_t(:, :, :)
      integer :: tinitx = 0
      logical :: diffuse = .false.
   end type marss_model

   type, public :: marss_constraint_block
      real(dp), allocatable :: fixed(:)
      real(dp), allocatable :: design(:, :)
      real(dp), allocatable :: start(:)
      character(len=marss_parameter_name_len), allocatable :: free_names(:)
   end type marss_constraint_block

   type, public :: marss_constraints
      type(marss_constraint_block) :: b
      type(marss_constraint_block) :: u
      type(marss_constraint_block) :: q
      type(marss_constraint_block) :: z
      type(marss_constraint_block) :: a
      type(marss_constraint_block) :: r
      type(marss_constraint_block) :: x0
      type(marss_constraint_block) :: v0
      type(marss_constraint_block) :: c
      type(marss_constraint_block) :: d
      type(marss_constraint_block) :: g
      type(marss_constraint_block) :: h
      type(marss_constraint_block) :: l
   end type marss_constraints

   type, public :: marss_kf_result
      logical :: ok = .false.
      integer :: info = 0
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp), allocatable :: x_pred(:, :)
      real(dp), allocatable :: p_pred(:, :, :)
      real(dp), allocatable :: p_pred_inf(:, :, :)
      real(dp), allocatable :: x_filt(:, :)
      real(dp), allocatable :: p_filt(:, :, :)
      real(dp), allocatable :: x_smooth(:, :)
      real(dp), allocatable :: p_smooth(:, :, :)
      real(dp), allocatable :: p_lag(:, :, :)
      real(dp), allocatable :: x0_smooth(:)
      real(dp), allocatable :: v0_smooth(:, :)
      real(dp), allocatable :: gain(:, :, :)
      real(dp), allocatable :: gain_inf(:, :, :)
      real(dp), allocatable :: innov(:, :)
      real(dp), allocatable :: sigma(:, :, :)
      real(dp), allocatable :: sigma_inf(:, :, :)
      integer :: diffuse_end = 0
      integer :: diffuse_component = 0
      integer :: remaining_diffuse_rank = 0
      logical :: observation_covariance_transformed = .false.
      logical :: lag_one_available = .true.
   end type marss_kf_result

   type, public :: marss_fit_result
      type(marss_model) :: model
      type(marss_kf_result) :: kf
      logical :: converged = .false.
      integer :: iterations = 0
      integer :: info = 0
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp), allocatable :: free_parameters(:)
   end type marss_fit_result

   type, public :: marss_hessian_result
      integer :: info = 0
      integer :: covariance_info = 0
      logical :: covariance_available = .false.
      logical :: had_nan_information = .false.
      character(len=16) :: method = ""
      real(dp), allocatable :: par_mean(:)
      real(dp), allocatable :: hessian(:, :)
      real(dp), allocatable :: par_sigma(:, :)
      character(len=marss_parameter_name_len), allocatable :: parameter_names(:)
   end type marss_hessian_result

   type, public :: marss_simulation
      real(dp), allocatable :: states(:, :, :)
      real(dp), allocatable :: data(:, :, :)
   end type marss_simulation

   type, public :: marss_cv_result
      real(dp), allocatable :: prediction(:, :)
      real(dp), allocatable :: se(:, :)
      integer, allocatable :: fold_id(:, :)
   end type marss_cv_result

   type, public :: marss_innov_boot_result
      real(dp), allocatable :: states(:, :, :)
      real(dp), allocatable :: data(:, :, :)
   end type marss_innov_boot_result

   type, public :: marss_residual_result
      integer :: info = 0
      real(dp), allocatable :: innovations(:, :)
      real(dp), allocatable :: standardized_innovations(:, :)
      real(dp), allocatable :: smoothed(:, :)
      real(dp), allocatable :: model_residuals(:, :)
      real(dp), allocatable :: state_residuals(:, :)
      real(dp), allocatable :: residuals(:, :)
      real(dp), allocatable :: var_residuals(:, :, :)
      real(dp), allocatable :: std_residuals(:, :)
      real(dp), allocatable :: marginal_residuals(:, :)
      real(dp), allocatable :: block_cholesky_residuals(:, :)
      real(dp), allocatable :: expected_observed_residuals(:, :)
      real(dp), allocatable :: var_observed_residuals(:, :, :)
   end type marss_residual_result

   type, public :: marss_hatyt_result
      logical :: ok = .false.
      integer :: info = 0
      real(dp), allocatable :: yt(:, :)
      real(dp), allocatable :: ot(:, :, :)
      real(dp), allocatable :: yxt(:, :, :)
      real(dp), allocatable :: yxtp(:, :, :)
      real(dp), allocatable :: var_yt(:, :, :)
      real(dp), allocatable :: var_expected_yt(:, :, :)
      real(dp), allocatable :: yxt_prev_smooth(:, :, :)
      real(dp), allocatable :: ytt1(:, :)
      real(dp), allocatable :: ott1(:, :, :)
      real(dp), allocatable :: var_ytt1(:, :, :)
      real(dp), allocatable :: var_expected_ytt1(:, :, :)
      real(dp), allocatable :: yxtt1(:, :, :)
      real(dp), allocatable :: ytt(:, :)
      real(dp), allocatable :: ott(:, :, :)
      real(dp), allocatable :: yxtt(:, :, :)
   end type marss_hatyt_result

end module marss_types
