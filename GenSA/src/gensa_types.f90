! SPDX-License-Identifier: GPL-2.0-only
module gensa_types
   use gensa_kinds, only : dp, i8
   implicit none
   private

   integer, parameter, public :: gensa_success = 0
   integer, parameter, public :: gensa_max_iterations = 1
   integer, parameter, public :: gensa_max_calls = 2
   integer, parameter, public :: gensa_max_time = 3
   integer, parameter, public :: gensa_threshold_reached = 4
   integer, parameter, public :: gensa_no_improvement = 5
   integer, parameter, public :: gensa_invalid_input = -1
   integer, parameter, public :: gensa_no_feasible_start = -2

   abstract interface
      function gensa_objective(x) result(value)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp) :: value
      end function gensa_objective

      logical function gensa_constraint(x)
         import dp
         real(dp), intent(in) :: x(:)
      end function gensa_constraint
   end interface

   type, public :: gensa_control
      integer :: maxit = 5000
      logical :: has_threshold = .false.
      real(dp) :: threshold_stop = 0.0_dp
      integer :: no_improvement_stop = -1
      logical :: smooth = .true.
      integer :: max_call = 10000000
      real(dp) :: max_time = huge(1.0_dp)
      real(dp) :: temperature = 5230.0_dp
      real(dp) :: visiting_param = 2.62_dp
      real(dp) :: acceptance_param = -5.0_dp
      logical :: verbose = .false.
      logical :: simple_function = .false.
      logical :: trace = .true.
      integer(i8) :: seed = -100377_i8
      integer :: markov_length = 0
      real(dp) :: temp_restart = 0.1_dp
      integer :: report = 100
      logical :: local_search = .true.
      integer :: local_maxit = 0
      real(dp) :: local_tolerance = 1.0e-7_dp
      integer :: max_constraint_attempts = 10000
   end type gensa_control

   type, public :: gensa_trace
      integer :: n = 0
      integer, allocatable :: step(:)
      real(dp), allocatable :: temperature(:)
      real(dp), allocatable :: current_value(:)
      real(dp), allocatable :: best_value(:)
   end type gensa_trace

   type, public :: gensa_result
      real(dp), allocatable :: par(:)
      real(dp) :: value = huge(1.0_dp)
      integer :: counts = 0
      integer :: iterations = 0
      integer :: status = gensa_invalid_input
      character(len=:), allocatable :: message
      type(gensa_trace) :: trace
   end type gensa_result

   public :: gensa_objective, gensa_constraint

end module gensa_types
