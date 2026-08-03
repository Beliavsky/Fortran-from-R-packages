! SPDX-License-Identifier: MIT
module fixedincome_types
   use fixedincome_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: FI_OK = 0
   integer, parameter, public :: FI_INVALID_ARGUMENT = 1
   integer, parameter, public :: FI_SIZE_MISMATCH = 2
   integer, parameter, public :: FI_OUT_OF_RANGE = 3
   integer, parameter, public :: FI_NOT_CONFIGURED = 4
   integer, parameter, public :: FI_NO_CONVERGENCE = 5
   integer, parameter, public :: FI_UNSUPPORTED_CALENDAR = 6

   integer, parameter, public :: COMPOUND_SIMPLE = 1
   integer, parameter, public :: COMPOUND_DISCRETE = 2
   integer, parameter, public :: COMPOUND_CONTINUOUS = 3

   integer, parameter, public :: UNIT_DAY = 1
   integer, parameter, public :: UNIT_MONTH = 2
   integer, parameter, public :: UNIT_YEAR = 3

   integer, parameter, public :: INTERP_NONE = 0
   integer, parameter, public :: INTERP_FLAT_FORWARD = 1
   integer, parameter, public :: INTERPOLATION_LINEAR = 2
   integer, parameter, public :: INTERP_LOG_LINEAR = 3
   integer, parameter, public :: INTERP_NATURAL_SPLINE = 4
   integer, parameter, public :: INTERP_HERMITE_SPLINE = 5
   integer, parameter, public :: INTERP_MONOTONE_SPLINE = 6
   integer, parameter, public :: INTERP_NELSON_SIEGEL = 7
   integer, parameter, public :: INTERP_NELSON_SIEGEL_SVENSSON = 8

   type, public :: term_t
      real(dp), allocatable :: value(:)
      integer, allocatable :: unit(:)
   contains
      procedure :: size => term_size
   end type term_t

   type, public :: daycount_t
      character(len=32) :: specification = 'actual/365'
      integer :: days_in_base = 365
   end type daycount_t

   type, public :: interpolation_t
      integer :: method = INTERP_NONE
      real(dp) :: parameters(6) = 0.0_dp
      logical :: propagate = .false.
   end type interpolation_t

   type, public :: spot_rate_t
      real(dp), allocatable :: rate(:)
      integer :: compounding = COMPOUND_DISCRETE
      type(daycount_t) :: daycount
      character(len=64) :: calendar = 'actual'
   contains
      procedure :: size => spot_rate_size
   end type spot_rate_t

   type, public :: spot_rate_curve_t
      real(dp), allocatable :: rate(:)
      real(dp), allocatable :: term_days(:)
      integer :: compounding = COMPOUND_DISCRETE
      type(daycount_t) :: daycount
      character(len=64) :: calendar = 'actual'
      integer :: reference_date = 0
      type(interpolation_t) :: interpolation
   contains
      procedure :: size => curve_size
   end type spot_rate_curve_t

   type, public :: forward_rate_t
      real(dp), allocatable :: rate(:)
      real(dp), allocatable :: interval_days(:)
      integer :: compounding = COMPOUND_DISCRETE
      type(daycount_t) :: daycount
      character(len=64) :: calendar = 'actual'
      integer :: reference_date = 0
   contains
      procedure :: size => forward_size
   end type forward_rate_t

   type, public :: fit_result_t
      type(interpolation_t) :: model
      real(dp) :: objective = huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = FI_NOT_CONFIGURED
   end type fit_result_t

contains

   pure integer function term_size(self) result(n)
      class(term_t), intent(in) :: self
      if (allocated(self%value)) then
         n = size(self%value)
      else
         n = 0
      end if
   end function term_size

   pure integer function spot_rate_size(self) result(n)
      class(spot_rate_t), intent(in) :: self
      if (allocated(self%rate)) then
         n = size(self%rate)
      else
         n = 0
      end if
   end function spot_rate_size

   pure integer function curve_size(self) result(n)
      class(spot_rate_curve_t), intent(in) :: self
      if (allocated(self%rate)) then
         n = size(self%rate)
      else
         n = 0
      end if
   end function curve_size

   pure integer function forward_size(self) result(n)
      class(forward_rate_t), intent(in) :: self
      if (allocated(self%rate)) then
         n = size(self%rate)
      else
         n = 0
      end if
   end function forward_size

end module fixedincome_types
