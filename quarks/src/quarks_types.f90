module quarks_types
   use quarks_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: quarks_ok = 0
   integer, parameter, public :: quarks_invalid_input = 1
   integer, parameter, public :: quarks_empty_tail = 2
   integer, parameter, public :: quarks_fit_failed = 3
   integer, parameter, public :: quarks_no_violations = 4

   integer, parameter, public :: method_plain = 1
   integer, parameter, public :: method_age = 2
   integer, parameter, public :: method_vwhs = 3
   integer, parameter, public :: method_fhs = 4

   integer, parameter, public :: volatility_ewma = 1
   integer, parameter, public :: volatility_garch = 2

   integer, parameter, public :: smooth_none = 0
   integer, parameter, public :: smooth_lpr = 1
   integer, parameter, public :: smooth_auto = 2

   type, public :: risk_result
      real(dp) :: var = 0.0_dp
      real(dp) :: es = 0.0_dp
      real(dp) :: p = 0.975_dp
      integer :: status = quarks_ok
      character(len=160) :: message = 'ok'
      integer :: volatility_model = volatility_ewma
      logical :: used_fallback = .false.
   end type risk_result

   type, public :: rollcast_result
      real(dp), allocatable :: var(:)
      real(dp), allocatable :: es(:)
      real(dp), allocatable :: xout(:)
      real(dp) :: p = 0.975_dp
      integer :: method = method_plain
      integer :: volatility_model = volatility_ewma
      integer :: smoothing = smooth_none
      integer :: nout = 0
      integer :: nwin = 0
      integer :: nboot = 0
      integer :: status = quarks_ok
      character(len=160) :: message = 'ok'
   end type rollcast_result

   type, public :: coverage_result
      real(dp) :: p = 0.975_dp
      real(dp) :: p_uc = 0.0_dp
      real(dp) :: p_ind = 0.0_dp
      real(dp) :: p_cc = 0.0_dp
      real(dp) :: lr_uc = 0.0_dp
      real(dp) :: lr_ind = 0.0_dp
      real(dp) :: lr_cc = 0.0_dp
      real(dp) :: confidence_level = 0.95_dp
      integer :: violations = 0
      integer :: n00 = 0
      integer :: n01 = 0
      integer :: n10 = 0
      integer :: n11 = 0
      integer :: status = quarks_ok
      character(len=160) :: message = 'ok'
   end type coverage_result

   type, public :: traffic_result
      real(dp) :: cumulative_probability = 0.0_dp
      real(dp) :: breach_probability = 0.025_dp
      integer :: violations = 0
      integer :: observations = 0
      integer :: status = quarks_ok
      character(len=160) :: message = 'ok'
   end type traffic_result

   type, public :: loss_result
      real(dp) :: lossfun1 = 0.0_dp
      real(dp) :: lossfun2 = 0.0_dp
      real(dp) :: lossfun3 = 0.0_dp
      real(dp) :: lossfun4 = 0.0_dp
      integer :: status = quarks_ok
      character(len=160) :: message = 'ok'
   end type loss_result

   type, public :: pl_result
      real(dp), allocatable :: pl(:)
      real(dp), allocatable :: weights(:,:)
      integer :: status = quarks_ok
      character(len=160) :: message = 'ok'
   end type pl_result

end module quarks_types
