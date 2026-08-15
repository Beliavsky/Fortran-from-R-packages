! Computational result/control types for the GAMLSS port.
! SPDX-License-Identifier: GPL-3.0-only
module gamlss_types
   use gamlss_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: GAMLSS_METHOD_RS = 1
   integer, parameter, public :: GAMLSS_METHOD_CG = 2
   integer, parameter, public :: GAMLSS_METHOD_MIXED = 3

   type, public :: gamlss_control_t
      real(dp) :: c_crit = 1.0e-3_dp
      integer :: n_cyc = 20
      integer :: inner_cyc = 50
      real(dp) :: inner_crit = 1.0e-5_dp
      real(dp) :: mu_step = 1.0_dp
      real(dp) :: sigma_step = 1.0_dp
      real(dp) :: nu_step = 1.0_dp
      real(dp) :: tau_step = 1.0_dp
      logical :: autostep = .true.
      integer :: mixed_rs_cycles = 1
      integer :: mixed_cg_cycles = 20
      ! Estimate the scalar quadratic-penalty multiplier during RS updates.
      ! This is used by random effects and ML-type penalized smoothers.
      logical :: estimate_lambda_mu = .false.
      logical :: estimate_lambda_sigma = .false.
      logical :: estimate_lambda_nu = .false.
      logical :: estimate_lambda_tau = .false.
      real(dp) :: lambda_min = 1.0e-7_dp
      real(dp) :: lambda_max = 1.0e7_dp
      real(dp) :: lambda_crit = 1.0e-6_dp
   end type gamlss_control_t

   type, public :: gamlss_parameter_result_t
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: eta(:)
      real(dp), allocatable :: covariance(:,:)
      real(dp) :: edf = 0.0_dp
      real(dp) :: penalty = 0.0_dp
      real(dp) :: lambda = 0.0_dp
   end type gamlss_parameter_result_t

   type, public :: gamlss_result_t
      type(gamlss_parameter_result_t) :: mu
      type(gamlss_parameter_result_t) :: sigma
      type(gamlss_parameter_result_t) :: nu
      type(gamlss_parameter_result_t) :: tau
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: case_deviance(:)
      real(dp) :: global_deviance = huge(1.0_dp)
      real(dp) :: penalized_deviance = huge(1.0_dp)
      real(dp) :: df_fit = 0.0_dp
      real(dp) :: df_residual = 0.0_dp
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: sbc = huge(1.0_dp)
      integer :: family = 0
      integer :: method = GAMLSS_METHOD_RS
      integer :: iterations = 0
      integer :: status = 0
      logical :: converged = .false.
   end type gamlss_result_t

end module gamlss_types
