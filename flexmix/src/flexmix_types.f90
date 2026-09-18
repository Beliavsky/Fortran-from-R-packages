! SPDX-License-Identifier: GPL-2.0-or-later
module flexmix_types
   use flexmix_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: flexmix_model_gaussian = 1
   integer, parameter, public :: flexmix_model_poisson = 2
   integer, parameter, public :: flexmix_model_binomial = 3
   integer, parameter, public :: flexmix_model_mvnorm = 4
   integer, parameter, public :: flexmix_model_mvbinary = 5
   integer, parameter, public :: flexmix_model_mvpois = 6
   integer, parameter, public :: flexmix_model_mvcombi = 7
   integer, parameter, public :: flexmix_model_lognormal = 8
   integer, parameter, public :: flexmix_model_exponential = 9
   integer, parameter, public :: flexmix_model_inverse_gaussian = 10
   integer, parameter, public :: flexmix_model_gamma = 11
   integer, parameter, public :: flexmix_model_weibull = 12
   integer, parameter, public :: flexmix_model_gamma_regression = 13
   integer, parameter, public :: flexmix_model_multinomial = 14
   integer, parameter, public :: flexmix_model_ziglm_poisson = 15
   integer, parameter, public :: flexmix_model_ziglm_binomial = 16
   integer, parameter, public :: flexmix_model_robust_gaussian = 17
   integer, parameter, public :: flexmix_model_robust_poisson = 18
   integer, parameter, public :: flexmix_model_factanal = 19
   integer, parameter, public :: flexmix_model_conditional_logit = 20
   integer, parameter, public :: flexmix_model_lmm = 21
   integer, parameter, public :: flexmix_model_lmer = 22
   integer, parameter, public :: flexmix_model_lmmc = 23
   integer, parameter, public :: flexmix_model_lmc = 24

   integer, parameter, public :: flexmix_class_weighted = 1
   integer, parameter, public :: flexmix_class_hard = 2
   integer, parameter, public :: flexmix_role_regular = 0
   integer, parameter, public :: flexmix_role_structural_zero = 1
   integer, parameter, public :: flexmix_role_robust_background = 2

   type, public :: flexmix_control
      integer :: iter_max = 200
      real(dp) :: minprior = 0.05_dp
      real(dp) :: tolerance = 1.0e-6_dp
      integer :: classify = flexmix_class_weighted
   end type flexmix_control

   type, public :: flexmix_result
      integer :: model_kind = 0
      integer :: n = 0
      integer :: p = 0
      integer :: d = 0
      integer :: nclass = 0
      integer :: factors = 0
      integer :: k = 0
      integer :: k0 = 0
      integer :: df = 0
      real(dp) :: effective_df = 0.0_dp
      integer :: iterations = 0
      integer :: status = 0
      logical :: converged = .false.
      logical :: diagonal_covariance = .true.
      logical :: truncated_binary = .false.
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp), allocatable :: prior(:)
      real(dp), allocatable :: posterior(:,:)
      real(dp), allocatable :: posterior_unscaled(:,:)
      real(dp), allocatable :: log_posterior_unscaled(:,:)
      integer, allocatable :: cluster(:)
      integer, allocatable :: size(:)
      integer, allocatable :: component_role(:)
      real(dp), allocatable :: beta(:,:)
      real(dp), allocatable :: multinomial_coef(:,:,:)
      real(dp), allocatable :: sigma(:)
      real(dp), allocatable :: center(:,:)
      real(dp), allocatable :: covariance(:,:,:)
      real(dp), allocatable :: factor_loadings(:,:,:)
      real(dp), allocatable :: uniqueness(:,:)
      real(dp), allocatable :: marginal_variance(:,:)
      real(dp), allocatable :: probability(:,:)
      real(dp), allocatable :: lambda(:,:)
      real(dp), allocatable :: shape(:)
      real(dp), allocatable :: rate(:)
      real(dp), allocatable :: scale(:)
      real(dp), allocatable :: concomitant_coef(:,:)
      real(dp), allocatable :: penalty_lambda(:)
      real(dp), allocatable :: smoothing_lambda(:)
      real(dp), allocatable :: component_effective_df(:)
      real(dp), allocatable :: random_covariance(:,:,:)
      real(dp), allocatable :: residual_variance(:)
      real(dp) :: penalty_alpha = 1.0_dp
      logical, allocatable :: binary_mask(:)
      logical, allocatable :: group_first(:)
      real(dp), allocatable :: case_weights(:)
      logical, allocatable :: parameter_design(:,:)
      integer, allocatable :: variance_group(:)
   end type flexmix_result

   type, public :: flexmix_step_result
      integer :: nrep = 0
      integer, allocatable :: k(:)
      real(dp), allocatable :: logliks(:,:)
      type(flexmix_result), allocatable :: models(:)
   end type flexmix_step_result

   type, public :: flexmix_boot_result
      integer :: nrep = 0
      integer, allocatable :: requested_k(:)
      real(dp), allocatable :: loglik(:,:)
      integer, allocatable :: fitted_k(:,:)
      logical, allocatable :: converged(:,:)
      type(flexmix_result), allocatable :: models(:,:)
   end type flexmix_boot_result

end module flexmix_types
