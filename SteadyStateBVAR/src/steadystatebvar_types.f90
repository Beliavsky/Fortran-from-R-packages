! SPDX-License-Identifier: GPL-3.0-or-later
! Derived from SteadyStateBVAR 0.2.0, Copyright (c) 2026 Mark Becker.
module steadystatebvar_types
   use r_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: ssbvar_success = 0
   integer, parameter, public :: ssbvar_invalid_input = 1
   integer, parameter, public :: ssbvar_not_ready = 2
   integer, parameter, public :: ssbvar_linalg_failure = 3
   integer, parameter, public :: ssbvar_unsupported = 4

   type, public :: bvar_setup
      integer :: n = 0
      integer :: k = 0
      integer :: p = 0
      integer :: q = 0
      integer :: n_free_params_a = 0
      character(len=24) :: deterministic = ""
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: x(:, :)
      real(dp), allocatable :: w(:, :)
      real(dp), allocatable :: q_lag(:, :)
      real(dp), allocatable :: d(:, :)
      real(dp), allocatable :: dt(:, :)
      real(dp), allocatable :: beta_ols(:, :)
      real(dp), allocatable :: sigma_u_ols(:, :)
      real(dp), allocatable :: psi_ols(:, :)
      real(dp), allocatable :: sigma_ar(:, :)
      real(dp), allocatable :: dummy(:)
   end type bvar_setup

   type, public :: bvar_sv_priors
      real(dp), allocatable :: theta_a(:)
      real(dp), allocatable :: omega_a(:, :)
      real(dp), allocatable :: theta_log_lambda_1(:)
      real(dp), allocatable :: omega_log_lambda_1(:, :)
      real(dp), allocatable :: alpha_phi(:)
      real(dp), allocatable :: beta_phi(:)
      real(dp), allocatable :: theta_gamma_0(:)
      real(dp), allocatable :: omega_gamma_0(:, :)
      real(dp), allocatable :: theta_gamma_1(:)
      real(dp), allocatable :: omega_gamma_1(:, :)
      integer :: m_phi = 0
      real(dp), allocatable :: v_phi(:, :)
   end type bvar_sv_priors

   type, public :: bvar_priors
      logical :: jeffreys = .true.
      logical :: sv = .false.
      character(len=8) :: sv_type = ""
      integer :: m = 0
      real(dp), allocatable :: theta_beta(:)
      real(dp), allocatable :: omega_beta(:, :)
      real(dp), allocatable :: theta_psi(:)
      real(dp), allocatable :: omega_psi(:, :)
      real(dp), allocatable :: sigma_ar(:, :)
      real(dp), allocatable :: v(:, :)
      real(dp), allocatable :: restriction(:, :)
      type(bvar_sv_priors) :: sv_priors
   end type bvar_priors

   type, public :: bvar_fit
      integer :: n_draws = 0
      integer :: h = 0
      logical :: homoscedastic = .true.
      character(len=8) :: sv_type = ""
      real(dp), allocatable :: d_pred(:, :)
      real(dp), allocatable :: beta(:, :, :)
      real(dp), allocatable :: psi(:, :, :)
      real(dp), allocatable :: sigma_u(:, :, :)
      real(dp), allocatable :: sigma_u_time(:, :, :, :)
      real(dp), allocatable :: sigma_u_pred(:, :, :, :)
      real(dp), allocatable :: y_pred(:, :, :)
      real(dp), allocatable :: a(:, :, :)
      real(dp), allocatable :: phi(:, :)
      real(dp), allocatable :: gamma_0(:, :)
      real(dp), allocatable :: gamma_1(:, :)
      real(dp), allocatable :: phi_cov(:, :, :)
      real(dp), allocatable :: log_lambda(:, :, :)
      real(dp), allocatable :: beta_mean(:, :)
      real(dp), allocatable :: beta_median(:, :)
      real(dp), allocatable :: psi_mean(:, :)
      real(dp), allocatable :: psi_median(:, :)
      real(dp), allocatable :: sigma_mean(:, :)
      real(dp), allocatable :: sigma_median(:, :)
      real(dp), allocatable :: sigma_time_mean(:, :, :)
      real(dp), allocatable :: sigma_time_median(:, :, :)
      real(dp), allocatable :: a_mean(:, :)
      real(dp), allocatable :: a_median(:, :)
      real(dp), allocatable :: phi_mean(:)
      real(dp), allocatable :: phi_median(:)
      real(dp), allocatable :: gamma_0_mean(:)
      real(dp), allocatable :: gamma_0_median(:)
      real(dp), allocatable :: gamma_1_mean(:)
      real(dp), allocatable :: gamma_1_median(:)
      real(dp), allocatable :: phi_cov_mean(:, :)
      real(dp), allocatable :: phi_cov_median(:, :)
   end type bvar_fit

   type, public :: bvar_model
      real(dp), allocatable :: data(:, :)
      type(bvar_setup) :: setup
      type(bvar_priors) :: priors
      type(bvar_fit) :: fit
      logical :: has_setup = .false.
      logical :: has_priors = .false.
      logical :: has_fit = .false.
   end type bvar_model

   type, public :: forecast_result
      real(dp), allocatable :: forecast(:, :)
      real(dp), allocatable :: lower(:, :)
      real(dp), allocatable :: upper(:, :)
   end type forecast_result

   type, public :: irf_result
      real(dp), allocatable :: center(:, :, :)
      real(dp), allocatable :: lower(:, :, :)
      real(dp), allocatable :: upper(:, :, :)
   end type irf_result

   type, public :: bvar_summary_result
      real(dp), allocatable :: beta(:, :)
      real(dp), allocatable :: psi(:, :)
      real(dp), allocatable :: sigma_u(:, :)
      real(dp), allocatable :: a(:, :)
      real(dp), allocatable :: phi(:)
      real(dp), allocatable :: gamma_0(:)
      real(dp), allocatable :: gamma_1(:)
      real(dp), allocatable :: phi_cov(:, :)
   end type bvar_summary_result

end module steadystatebvar_types
