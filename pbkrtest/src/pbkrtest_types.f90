! SPDX-License-Identifier: GPL-2.0-or-later
module pbkrtest_types
   use r_kinds, only : dp
   use numderiv_callbacks, only : scalar_real_function, vector_real_function
   implicit none
   private

   integer, parameter, public :: pbkr_success = 0
   integer, parameter, public :: pbkr_invalid_shape = -1
   integer, parameter, public :: pbkr_invalid_argument = -2
   integer, parameter, public :: pbkr_linalg_failure = -3
   integer, parameter, public :: pbkr_numderiv_failure = -4

   type, public :: random_sigma_term_t
      integer :: n_levels = 0
      real(dp), allocatable :: z(:, :)
      real(dp), allocatable :: covariance(:, :)
   end type random_sigma_term_t

   type, public :: sigma_g_result_t
      real(dp), allocatable :: sigma(:, :)
      real(dp), allocatable :: g(:, :, :)
      real(dp), allocatable :: gamma(:)
   end type sigma_g_result_t

   type, public :: vcov_adjustment_t
      real(dp), allocatable :: phi_adjusted(:, :)
      real(dp), allocatable :: p_matrices(:, :, :)
      real(dp), allocatable :: w(:, :)
      real(dp), allocatable :: information(:, :)
      real(dp) :: condition = 0.0_dp
   end type vcov_adjustment_t

   type, public :: kr_result_t
      integer :: ndf = 0
      real(dp) :: ddf = 0.0_dp
      real(dp) :: f_stat = 0.0_dp
      real(dp) :: p_value = 0.0_dp
      real(dp) :: f_scaling = 0.0_dp
      real(dp) :: f_stat_unscaled = 0.0_dp
      real(dp) :: p_value_unscaled = 0.0_dp
      real(dp) :: wald_unadjusted = 0.0_dp
      real(dp) :: a1 = 0.0_dp
      real(dp) :: a2 = 0.0_dp
      real(dp) :: v0 = 0.0_dp
      real(dp) :: v1 = 0.0_dp
      real(dp) :: v2 = 0.0_dp
      real(dp) :: rho = 0.0_dp
   end type kr_result_t

   type, public :: satterthwaite_result_t
      integer :: ndf = 0
      real(dp) :: ddf = 0.0_dp
      real(dp) :: f_stat = 0.0_dp
      real(dp) :: p_value = 0.0_dp
      real(dp), allocatable :: nu(:)
   end type satterthwaite_result_t

   type, public :: auxiliary_callbacks_t
      procedure(scalar_real_function), pointer, nopass :: deviance => null()
      procedure(vector_real_function), pointer, nopass :: covbeta_vector => null()
   end type auxiliary_callbacks_t

   type, public :: auxiliary_result_t
      real(dp), allocatable :: vcov_varpar(:, :)
      real(dp), allocatable :: jacobian(:, :, :)
      real(dp), allocatable :: hessian(:, :)
      real(dp), allocatable :: hessian_eigenvalues(:)
      integer :: negative_eigenvalues = 0
      integer :: near_zero_eigenvalues = 0
   end type auxiliary_result_t

   type, public :: lrt_result_t
      real(dp) :: statistic = 0.0_dp
      integer :: df = 0
      real(dp) :: p_value = 0.0_dp
   end type lrt_result_t

   type, public :: bootstrap_result_t
      integer :: nsim = 0
      integer :: npos = 0
      integer :: n_extreme = 0
      real(dp) :: mean_positive = 0.0_dp
      real(dp) :: variance_positive = 0.0_dp
      real(dp) :: p_chisq = 0.0_dp
      real(dp) :: p_bootstrap = 0.0_dp
      real(dp) :: p_bootstrap_all = 0.0_dp
      real(dp) :: standard_error = 0.0_dp
      real(dp) :: ci_low = 0.0_dp
      real(dp) :: ci_high = 0.0_dp
      real(dp) :: bartlett_statistic = 0.0_dp
      real(dp) :: p_bartlett = 0.0_dp
      real(dp) :: gamma_scale = 0.0_dp
      real(dp) :: gamma_shape = 0.0_dp
      real(dp) :: p_gamma = 0.0_dp
      real(dp) :: f_statistic = 0.0_dp
      real(dp) :: f_ddf = 0.0_dp
      real(dp) :: p_f = 0.0_dp
   end type bootstrap_result_t

   public :: dp

end module pbkrtest_types
