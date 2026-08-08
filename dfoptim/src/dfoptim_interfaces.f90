! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of computational code from the R package dfoptim.
! Original authors: Ravi Varadhan, Hans W. Borchers, and Vincent Bechard.

module dfoptim_interfaces
   use dfoptim_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: dfoptim_success = 0
   integer, parameter, public :: dfoptim_max_evaluations = 1
   integer, parameter, public :: dfoptim_stagnation = 2
   integer, parameter, public :: dfoptim_invalid_input = -1
   integer, parameter, public :: dfoptim_nonfinite_objective = -2
   integer, parameter, public :: dfoptim_cancelled = -3

   integer, parameter, public :: mads_poll_lite = 1
   integer, parameter, public :: mads_poll_full = 2

   abstract interface
      function dfoptim_objective(x, user_data) result(value)
         import :: dp
         real(dp), intent(in) :: x(:)
         class(*), intent(inout), optional :: user_data
         real(dp) :: value
      end function dfoptim_objective

      subroutine dfoptim_monitor(x, value, iteration, evaluations, stop, user_data)
         import :: dp
         real(dp), intent(in) :: x(:)
         real(dp), intent(in) :: value
         integer, intent(in) :: iteration
         integer, intent(in) :: evaluations
         logical, intent(out) :: stop
         class(*), intent(inout), optional :: user_data
      end subroutine dfoptim_monitor
   end interface

   type, public :: hj_control_t
      real(dp) :: tol = 1.0e-6_dp
      integer :: maxfeval = 500000000
      logical :: maximize = .false.
      real(dp) :: target = huge(1.0_dp)
      logical :: trace = .false.
      integer :: seed = 1138
   end type hj_control_t

   type, public :: nmk_control_t
      real(dp) :: tol = 1.0e-6_dp
      integer :: maxfeval = 0
      logical :: regular_simplex = .true.
      logical :: maximize = .false.
      integer :: max_restarts = 3
      logical :: trace = .false.
   end type nmk_control_t

   type, public :: mads_control_t
      logical :: trace = .false.
      real(dp) :: tol = 1.0e-6_dp
      integer :: maxfeval = 10000
      logical :: maximize = .false.
      integer :: poll_style = mads_poll_lite
      real(dp) :: delta_init = 0.01_dp
      real(dp) :: expand = 4.0_dp
      integer :: line_search = 20
      integer :: seed = 1138
   end type mads_control_t

   type, public :: mads_log_t
      integer, allocatable :: evaluations(:)
      real(dp), allocatable :: delta(:)
      integer, allocatable :: search_success(:)
      real(dp), allocatable :: value(:)
   end type mads_log_t

   type, public :: dfoptim_result_t
      real(dp), allocatable :: x(:)
      real(dp) :: value = huge(1.0_dp)
      integer :: convergence = dfoptim_invalid_input
      integer :: feval = 0
      integer :: niter = 0
      integer :: restarts = 0
      real(dp) :: final_mesh = huge(1.0_dp)
      character(len=:), allocatable :: message
      type(mads_log_t) :: log
   end type dfoptim_result_t

   public :: dfoptim_objective, dfoptim_monitor

end module dfoptim_interfaces
