module lme4_types
   use lme4_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: family_binomial = 1
   integer, parameter, public :: family_poisson = 2
   integer, parameter, public :: family_gamma = 3
   integer, parameter, public :: family_inverse_gaussian = 4
   integer, parameter, public :: family_negative_binomial = 5

   integer, parameter, public :: covariance_unstructured = 1
   integer, parameter, public :: covariance_diagonal = 2
   integer, parameter, public :: covariance_compound_symmetry = 3
   integer, parameter, public :: covariance_ar1 = 4

   type, public :: random_term_t
      real(dp), allocatable :: z(:,:)
      integer, allocatable :: group(:)
      integer :: n_levels = 0
      character(len=:), allocatable :: name
      integer :: covariance_structure = covariance_unstructured
   contains
      procedure :: validate => validate_random_term
      procedure :: n_coefficients => term_n_coefficients
      procedure :: n_random_effects => term_n_random_effects
      procedure :: n_parameters => term_n_parameters
   end type random_term_t

   type, public :: covariance_block_t
      character(len=:), allocatable :: name
      real(dp), allocatable :: covariance(:,:)
      real(dp), allocatable :: sdcor(:,:)
      integer :: n_levels = 0
   end type covariance_block_t

   type, public :: lmm_control_t
      integer :: maxfun = 5000
      real(dp) :: tolerance = 1.0e-7_dp
      real(dp) :: lower_log_sd = -8.0_dp
      real(dp) :: upper_log_sd = 4.0_dp
      real(dp) :: lower_offdiag = -5.0_dp
      real(dp) :: upper_offdiag = 5.0_dp
      logical :: verbose = .false.
   end type lmm_control_t

   type, public :: glmm_control_t
      integer :: maxfun = 5000
      integer :: max_pirls = 100
      real(dp) :: tolerance = 1.0e-7_dp
      real(dp) :: pirls_tolerance = 1.0e-8_dp
      real(dp) :: lower_log_sd = -8.0_dp
      real(dp) :: upper_log_sd = 3.0_dp
      real(dp) :: lower_offdiag = -5.0_dp
      real(dp) :: upper_offdiag = 5.0_dp
      real(dp) :: lower_beta = -30.0_dp
      real(dp) :: upper_beta = 30.0_dp
      real(dp) :: lower_log_dispersion = -4.605170185988091_dp
      real(dp) :: upper_log_dispersion = 9.210340371976184_dp
      integer :: max_profile = 40
      logical :: verbose = .false.
   end type glmm_control_t

   type, public :: lmm_result_t
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: u(:)
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: vcov_beta(:,:)
      type(covariance_block_t), allocatable :: varcorr(:)
      integer, allocatable :: term_offsets(:)
      real(dp) :: sigma = 0.0_dp
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: deviance = huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: evaluations = 0
      integer :: status = 0
      logical :: converged = .false.
      logical :: reml = .true.
      character(len=:), allocatable :: message
   end type lmm_result_t

   type, public :: glmm_result_t
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: u(:)
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: linear_predictor(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: vcov_beta(:,:)
      type(covariance_block_t), allocatable :: varcorr(:)
      integer, allocatable :: term_offsets(:)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: deviance = huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      real(dp) :: dispersion = 1.0_dp
      integer :: family = 0
      integer :: quadrature_order = 1
      integer :: evaluations = 0
      integer :: pirls_iterations = 0
      integer :: status = 0
      logical :: converged = .false.
      character(len=:), allocatable :: message
   end type glmm_result_t


   type, public :: nlmm_control_t
      integer :: maxfun = 5000
      integer :: max_mode_iterations = 60
      real(dp) :: tolerance = 1.0e-7_dp
      real(dp) :: mode_tolerance = 1.0e-8_dp
      real(dp) :: lower_beta = -30.0_dp
      real(dp) :: upper_beta = 30.0_dp
      real(dp) :: lower_log_sd = -8.0_dp
      real(dp) :: upper_log_sd = 4.0_dp
      real(dp) :: lower_offdiag = -5.0_dp
      real(dp) :: upper_offdiag = 5.0_dp
      real(dp) :: lower_log_sigma = -10.0_dp
      real(dp) :: upper_log_sigma = 5.0_dp
      logical :: verbose = .false.
   end type nlmm_control_t

   type, public :: nlmm_result_t
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: u(:,:)
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: vcov_beta(:,:)
      real(dp), allocatable :: covariance(:,:)
      real(dp) :: sigma = 0.0_dp
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: deviance = huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: evaluations = 0
      integer :: status = 0
      logical :: converged = .false.
      character(len=:), allocatable :: message
   end type nlmm_result_t

   type, public :: bootstrap_result_t
      real(dp), allocatable :: beta(:,:)
      real(dp), allocatable :: theta(:,:)
      real(dp), allocatable :: scale(:)
      real(dp), allocatable :: log_likelihood(:)
      logical, allocatable :: converged(:)
      integer :: successful = 0
   end type bootstrap_result_t

   type, public :: profile_result_t
      real(dp), allocatable :: estimate(:)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: upper(:)
      real(dp), allocatable :: standard_error(:)
      real(dp) :: level = 0.95_dp
      integer :: status = 0
      character(len=:), allocatable :: message
   end type profile_result_t

   type, public :: lm_list_result_t
      real(dp), allocatable :: coefficients(:,:)
      real(dp), allocatable :: sigma(:)
      real(dp), allocatable :: log_likelihood(:)
      integer, allocatable :: observations(:)
      logical, allocatable :: converged(:)
   end type lm_list_result_t

   type, public :: influence_result_t
      real(dp), allocatable :: dfbeta(:,:)
      real(dp), allocatable :: cooks_distance(:)
      real(dp), allocatable :: deleted_log_likelihood(:)
      logical, allocatable :: converged(:)
   end type influence_result_t

   type, public :: gh_rule_t
      real(dp), allocatable :: nodes(:)
      real(dp), allocatable :: weights(:)
      real(dp), allocatable :: log_density(:)
   end type gh_rule_t

contains

   subroutine validate_random_term(self, n, ok, message)
      class(random_term_t), intent(inout) :: self
      integer, intent(in) :: n
      logical, intent(out) :: ok
      character(len=:), allocatable, intent(out) :: message

      ok = .false.
      if (.not. allocated(self%z) .or. .not. allocated(self%group)) then
         message = 'random term requires allocated z and group arrays'
         return
      end if
      if (size(self%z,1) /= n .or. size(self%group) /= n) then
         message = 'random term row counts must match the response'
         return
      end if
      if (size(self%z,2) < 1) then
         message = 'random term must have at least one coefficient'
         return
      end if
      if (any(self%group < 1)) then
         message = 'group indices must be positive integers'
         return
      end if
      if (self%n_levels <= 0) self%n_levels = maxval(self%group)
      if (maxval(self%group) > self%n_levels) then
         message = 'group index exceeds n_levels'
         return
      end if
      if (.not. allocated(self%name)) self%name = 'group'
      select case (self%covariance_structure)
      case (covariance_unstructured, covariance_diagonal)
      case (covariance_compound_symmetry, covariance_ar1)
         if (size(self%z,2) < 2) then
            message = 'compound-symmetry and AR(1) covariance require at least two coefficients'
            return
         end if
      case default
         message = 'unknown covariance structure'
         return
      end select
      ok = .true.
      message = 'ok'
   end subroutine validate_random_term

   integer function term_n_coefficients(self) result(q)
      class(random_term_t), intent(in) :: self
      if (allocated(self%z)) then
         q = size(self%z,2)
      else
         q = 0
      end if
   end function term_n_coefficients

   integer function term_n_random_effects(self) result(nr)
      class(random_term_t), intent(in) :: self
      nr = self%n_levels * self%n_coefficients()
   end function term_n_random_effects

   integer function term_n_parameters(self) result(nt)
      class(random_term_t), intent(in) :: self
      integer :: q
      q = self%n_coefficients()
      select case (self%covariance_structure)
      case (covariance_unstructured)
         nt = q * (q + 1) / 2
      case (covariance_diagonal)
         nt = q
      case (covariance_compound_symmetry, covariance_ar1)
         nt = 2
      case default
         nt = 0
      end select
   end function term_n_parameters

end module lme4_types
