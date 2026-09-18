module gamm4_types
   use gamm4_kinds, only : dp
   use mgcv, only : smooth_spec_t
   implicit none
   private

   integer, parameter, public :: gamm4_family_gaussian = 0

   type, public :: gamm4_smooth_t
      real(dp), allocatable :: basis(:, :)
      type(smooth_spec_t) :: spec
      character(len=:), allocatable :: label
   end type gamm4_smooth_t

   type, public :: gamm4_control_t
      integer :: max_outer = 12
      integer :: max_coordinate = 48
      integer :: max_pirls = 100
      real(dp) :: tolerance = 1.0e-5_dp
      real(dp) :: pirls_tolerance = 1.0e-8_dp
      real(dp) :: lower_log_sd = -8.0_dp
      real(dp) :: upper_log_sd = 4.0_dp
      real(dp) :: lower_offdiag = -5.0_dp
      real(dp) :: upper_offdiag = 5.0_dp
      real(dp) :: eigen_tolerance = sqrt(epsilon(1.0_dp))
      logical :: trace = .false.
   end type gamm4_control_t

   type, public :: gamm4_result_t
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: edf(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: linear_predictor(:)
      real(dp), allocatable :: conditional_fitted(:)
      real(dp), allocatable :: conditional_linear_predictor(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: sp(:)
      real(dp), allocatable :: smooth_sd(:)
      real(dp), allocatable :: variance_parameters(:)
      real(dp), allocatable :: random_effects(:)
      real(dp) :: scale = 1.0_dp
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: deviance = huge(1.0_dp)
      integer :: family = gamm4_family_gaussian
      integer :: iterations = 0
      integer :: status = 0
      logical :: converged = .false.
      logical :: reml = .true.
      character(len=:), allocatable :: method
      character(len=:), allocatable :: message
   end type gamm4_result_t

   type, public :: gamm4_vb_result_t
      real(dp), allocatable :: vb(:, :)
      real(dp), allocatable :: xvx(:, :)
      real(dp), allocatable :: r(:, :)
      integer :: status = 0
   end type gamm4_vb_result_t

end module gamm4_types
