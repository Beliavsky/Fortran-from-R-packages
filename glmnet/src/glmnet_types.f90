! SPDX-License-Identifier: GPL-2.0-only
module glmnet_types
   use glmnet_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: glmnet_family_gaussian = 1
   integer, parameter, public :: glmnet_family_binomial = 2
   integer, parameter, public :: glmnet_family_poisson = 3
   integer, parameter, public :: glmnet_family_multinomial = 4
   integer, parameter, public :: glmnet_family_mgaussian = 5
   integer, parameter, public :: glmnet_family_cox = 6

   type, public :: glmnet_control_type
      real(dp) :: alpha = 1.0_dp
      integer :: nlambda = 100
      real(dp) :: lambda_min_ratio = -1.0_dp
      real(dp) :: threshold = 1.0e-7_dp
      integer :: max_iterations = 100000
      integer :: max_outer_iterations = 100
      integer :: max_active = huge(1)
      real(dp) :: fractional_deviance = 1.0e-5_dp
      real(dp) :: deviance_max = 0.999_dp
      integer :: minimum_lambda_count = 5
      logical :: standardize = .true.
      logical :: intercept = .true.
      logical :: grouped = .false.
      logical :: trace = .false.
      real(dp) :: probability_min = 1.0e-9_dp
      real(dp) :: eta_max = 30.0_dp
      real(dp) :: step_min = 1.0e-12_dp
   end type glmnet_control_type

   type, public :: glmnet_path_result
      character(len=16) :: family = ''
      integer :: family_code = 0
      integer :: status = 0
      integer :: nobs = 0
      integer :: nvars = 0
      integer :: nout = 0
      integer :: nlambda = 0
      integer :: npasses = 0
      logical :: standardize = .false.
      logical :: intercept_fitted = .false.
      logical :: efron = .false.
      real(dp) :: nulldev = 0.0_dp
      real(dp), allocatable :: lambda(:)
      real(dp), allocatable :: intercept(:,:)
      real(dp), allocatable :: beta(:,:,:)
      real(dp), allocatable :: dev_ratio(:)
      real(dp), allocatable :: objective(:)
      integer, allocatable :: df(:)
      integer, allocatable :: iterations(:)
      logical, allocatable :: converged(:)
      real(dp), allocatable :: x_mean(:)
      real(dp), allocatable :: x_scale(:)
      real(dp), allocatable :: class_levels(:)
   end type glmnet_path_result

   type, public :: glmnet_cv_result
      integer :: status = 0
      integer :: nfolds = 0
      integer :: index_min = 0
      integer :: index_1se = 0
      real(dp) :: lambda_min = 0.0_dp
      real(dp) :: lambda_1se = 0.0_dp
      character(len=16) :: measure = ''
      type(glmnet_path_result) :: fit
      real(dp), allocatable :: cv_mean(:)
      real(dp), allocatable :: cv_sd(:)
      real(dp), allocatable :: fold_values(:,:)
      real(dp), allocatable :: predictions(:,:,:)
      integer, allocatable :: fold_id(:)
   end type glmnet_cv_result

   type, public :: glmnet_assessment_result
      integer :: status = 0
      character(len=16) :: measure = ''
      real(dp), allocatable :: value(:)
      real(dp), allocatable :: secondary(:)
   end type glmnet_assessment_result

   type, public :: glmnet_roc_result
      integer :: status = 0
      real(dp) :: auc = 0.0_dp
      real(dp), allocatable :: threshold(:)
      real(dp), allocatable :: false_positive_rate(:)
      real(dp), allocatable :: true_positive_rate(:)
   end type glmnet_roc_result

   type, public :: glmnet_survival_data
      real(dp), allocatable :: start(:)
      real(dp), allocatable :: stop(:)
      integer, allocatable :: event(:)
      integer, allocatable :: strata(:)
   end type glmnet_survival_data

   type, public :: glmnet_sparse_csc
      integer :: nrow = 0
      integer :: ncol = 0
      real(dp), allocatable :: values(:)
      integer, allocatable :: row_index(:)
      integer, allocatable :: col_pointer(:)
   end type glmnet_sparse_csc

   abstract interface
      subroutine glmnet_family_working_interface(y, eta, base_weight, working, &
         irls_weight, deviance, status)
         import dp
         real(dp), intent(in) :: y(:), eta(:), base_weight(:)
         real(dp), intent(out) :: working(:), irls_weight(:)
         real(dp), intent(out) :: deviance
         integer, intent(out) :: status
      end subroutine glmnet_family_working_interface
   end interface
   public :: glmnet_family_working_interface, family_name

contains
   pure function family_name(code) result(name)
      integer, intent(in) :: code
      character(len=16) :: name
      select case (code)
      case (glmnet_family_gaussian)
         name = 'gaussian'
      case (glmnet_family_binomial)
         name = 'binomial'
      case (glmnet_family_poisson)
         name = 'poisson'
      case (glmnet_family_multinomial)
         name = 'multinomial'
      case (glmnet_family_mgaussian)
         name = 'mgaussian'
      case (glmnet_family_cox)
         name = 'cox'
      case default
         name = 'unknown'
      end select
   end function family_name
end module glmnet_types
