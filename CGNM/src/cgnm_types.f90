! SPDX-License-Identifier: MIT
module cgnm_types
   use cgnm_kinds, only : dp
   implicit none
   private

   public :: cgnm_model, cgnm_options, cgnm_problem, cgnm_result
   public :: cgnm_init_problem

   abstract interface
      subroutine cgnm_model(x, y, ierr)
         import :: dp
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: y(:)
         integer, intent(out) :: ierr
      end subroutine cgnm_model
   end interface

   type :: cgnm_options
      integer :: num_minimizers = 250
      integer :: num_iterations = 25
      real(dp) :: initial_lambda = 1.0_dp
      real(dp) :: gamma = 2.0_dp
      logical :: gamma_auto = .false.
      logical :: stay_in_initial_range = .false.
      integer :: seed = 12345
      integer :: max_initial_draws = 10000
      integer :: kmeans_max_iter = 100
   end type cgnm_options

   type :: cgnm_problem
      procedure(cgnm_model), pointer, nopass :: model => null()
      integer :: npar = 0
      integer :: nobs_model = 0
      integer :: nobs_total = 0
      real(dp), allocatable :: target(:)
      real(dp), allocatable :: initial_lower(:), initial_upper(:)
      real(dp), allocatable :: lower_bound(:), upper_bound(:)
      logical, allocatable :: has_lower(:), has_upper(:)
      real(dp), allocatable :: mo_weights(:), mo_values(:)
      logical, allocatable :: keep_initial_distribution(:)
   end type cgnm_problem

   type :: cgnm_result
      real(dp), allocatable :: x(:,:)              ! internal coordinates
      real(dp), allocatable :: theta(:,:)          ! physical parameters
      real(dp), allocatable :: y(:,:)
      real(dp), allocatable :: initial_x(:,:)
      real(dp), allocatable :: initial_theta(:,:)
      real(dp), allocatable :: initial_y(:,:)
      real(dp), allocatable :: residual_history(:,:)
      real(dp), allocatable :: lambda_history(:,:)
      real(dp), allocatable :: target_matrix(:,:)
      real(dp), allocatable :: weight_matrix(:,:)
      integer :: iterations = 0
      integer :: status = 0
      character(len=160) :: message = ''
   end type cgnm_result

contains

   subroutine cgnm_init_problem(prob, model, target, initial_lower, initial_upper, &
                                lower_bound, upper_bound, lower_active, upper_active, &
                                mo_weights, mo_values, keep_initial_distribution, ierr)
      type(cgnm_problem), intent(out) :: prob
      procedure(cgnm_model) :: model
      real(dp), intent(in) :: target(:), initial_lower(:), initial_upper(:)
      real(dp), intent(in), optional :: lower_bound(:), upper_bound(:)
      logical, intent(in), optional :: lower_active(:), upper_active(:)
      real(dp), intent(in), optional :: mo_weights(:), mo_values(:)
      logical, intent(in), optional :: keep_initial_distribution(:)
      integer, intent(out), optional :: ierr
      integer :: n, nm

      if (present(ierr)) ierr = 0
      n = size(initial_lower)
      if (size(initial_upper) /= n .or. n < 1) then
         if (present(ierr)) ierr = 1
         return
      end if
      prob%model => model
      prob%npar = n
      prob%nobs_model = size(target)
      prob%target = target
      prob%initial_lower = initial_lower
      prob%initial_upper = initial_upper
      allocate(prob%lower_bound(n), prob%upper_bound(n), prob%has_lower(n), prob%has_upper(n))
      prob%lower_bound = 0.0_dp
      prob%upper_bound = 0.0_dp
      prob%has_lower = .false.
      prob%has_upper = .false.
      if (present(lower_bound)) then
         if (size(lower_bound) /= n) then
            if (present(ierr)) ierr = 2
            return
         end if
         prob%lower_bound = lower_bound
         prob%has_lower = .true.
         if (present(lower_active)) then
            if (size(lower_active) /= n) then
               if (present(ierr)) ierr = 7
               return
            end if
            prob%has_lower = lower_active
         end if
      end if
      if (present(upper_bound)) then
         if (size(upper_bound) /= n) then
            if (present(ierr)) ierr = 3
            return
         end if
         prob%upper_bound = upper_bound
         prob%has_upper = .true.
         if (present(upper_active)) then
            if (size(upper_active) /= n) then
               if (present(ierr)) ierr = 8
               return
            end if
            prob%has_upper = upper_active
         end if
      end if
      allocate(prob%mo_weights(n), prob%mo_values(n), prob%keep_initial_distribution(n))
      prob%mo_weights = 0.0_dp
      prob%mo_values = 0.5_dp * (initial_lower + initial_upper)
      prob%keep_initial_distribution = .false.
      if (present(mo_weights)) then
         if (size(mo_weights) /= n) then
            if (present(ierr)) ierr = 4
            return
         end if
         prob%mo_weights = mo_weights
      end if
      if (present(mo_values)) then
         if (size(mo_values) /= n) then
            if (present(ierr)) ierr = 5
            return
         end if
         prob%mo_values = mo_values
      end if
      if (present(keep_initial_distribution)) then
         if (size(keep_initial_distribution) /= n) then
            if (present(ierr)) ierr = 6
            return
         end if
         prob%keep_initial_distribution = keep_initial_distribution
      end if
      nm = count(abs(prob%mo_weights) > 0.0_dp)
      prob%nobs_total = prob%nobs_model + nm
   end subroutine cgnm_init_problem

end module cgnm_types
