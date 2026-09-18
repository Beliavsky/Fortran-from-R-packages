module grpreg_types
   use grpreg_kinds, only : dp
   implicit none
   private

   type, public :: grpreg_fit_type
      character(len=16) :: family = 'gaussian'
      character(len=16) :: penalty = 'grLasso'
      integer :: n = 0
      integer :: p = 0
      integer :: ngroups = 0
      integer :: nlambda = 0
      real(dp) :: alpha = 1.0_dp
      real(dp) :: gamma = 3.0_dp
      real(dp) :: tau = 1.0_dp/3.0_dp
      real(dp), allocatable :: lambda(:)
      real(dp), allocatable :: intercept(:)
      real(dp), allocatable :: beta(:,:)
      real(dp), allocatable :: deviance(:)
      real(dp), allocatable :: df(:)
      integer, allocatable :: iter(:)
      real(dp), allocatable :: eta(:,:)
      integer, allocatable :: group(:)
      real(dp), allocatable :: group_multiplier(:)
      real(dp), allocatable :: x_center(:)
      real(dp), allocatable :: x_scale(:)
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: time(:)
      real(dp), allocatable :: fail(:)
   end type grpreg_fit_type

   type, public :: grpreg_cv_type
      type(grpreg_fit_type) :: fit
      integer :: nfolds = 0
      integer :: min_index = 0
      real(dp) :: lambda_min = 0.0_dp
      real(dp) :: null_deviance = 0.0_dp
      real(dp), allocatable :: lambda(:)
      real(dp), allocatable :: cve(:)
      real(dp), allocatable :: cvse(:)
      real(dp), allocatable :: prediction_error(:)
      integer, allocatable :: fold(:)
      real(dp), allocatable :: cv_prediction(:,:)
   end type grpreg_cv_type

   type, public :: grpreg_selection_type
      integer :: index = 0
      real(dp) :: lambda = 0.0_dp
      real(dp) :: df = 0.0_dp
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: criterion(:)
   end type grpreg_selection_type

   type, public :: mfdr_result_type
      real(dp), allocatable :: ef(:)
      integer, allocatable :: selected(:)
      real(dp), allocatable :: mfdr(:)
   end type mfdr_result_type

   type, public :: spline_expansion_type
      real(dp), allocatable :: x(:,:)
      integer, allocatable :: group(:)
      real(dp), allocatable :: knots(:,:)
      real(dp), allocatable :: boundary(:,:)
      integer :: degree = 3
      integer :: df = 3
      character(len=2) :: spline_type = 'bs'
   end type spline_expansion_type

end module grpreg_types
