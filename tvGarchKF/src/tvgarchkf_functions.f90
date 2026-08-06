! Modern Fortran translation of the computational core of tvGarchKF.
! Copyright (C) 2026 OpenAI
! SPDX-License-Identifier: GPL-3.0-or-later
module tvgarchkf_functions
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use fgarch_kinds, only : dp
   use tvgarchkf_types
   implicit none
   private

   public :: evaluate_tv_function, polynomial_values, nonlinear_values
   public :: trigonometric_values, identity_argument, three_one_minus_log_argument

contains

   function polynomial_values(n, coefficients) result(values)
      integer, intent(in) :: n
      real(dp), intent(in) :: coefficients(:)
      real(dp), allocatable :: values(:)
      real(dp) :: u
      integer :: i, j

      allocate(values(max(0,n)))
      values = 0.0_dp
      if (n <= 0) return
      do i = 1, n
         u = real(i,dp)/real(n,dp)
         do j = 1, size(coefficients)
            values(i) = values(i)+coefficients(j)*u**(j-1)
         end do
      end do
   end function polynomial_values

   function nonlinear_values(n, coefficients, exponents) result(values)
      integer, intent(in) :: n
      real(dp), intent(in) :: coefficients(:), exponents(:)
      real(dp), allocatable :: values(:)
      real(dp) :: u
      integer :: i, j, m

      allocate(values(max(0,n)))
      values = 0.0_dp
      if (n <= 0) return
      m = min(size(coefficients),size(exponents))
      do i = 1, n
         u = real(i,dp)/real(n,dp)
         do j = 1, m
            values(i) = values(i)+coefficients(j)*u**exponents(j)
         end do
      end do
   end function nonlinear_values

   function trigonometric_values(n, coefficients, trig_kind, argument_kind, custom_argument) result(values)
      integer, intent(in) :: n, trig_kind, argument_kind
      real(dp), intent(in) :: coefficients(:)
      procedure(tv_argument_function), pointer, intent(in), optional :: custom_argument
      real(dp), allocatable :: values(:)
      real(dp) :: u, argument, basis
      integer :: i

      allocate(values(max(0,n)))
      values = 0.0_dp
      if (n <= 0 .or. size(coefficients) == 0) return
      do i = 1, n
         u = real(i,dp)/real(n,dp)
         select case (argument_kind)
         case (arg_identity)
            argument = identity_argument(u)
         case (arg_three_one_minus_log)
            argument = three_one_minus_log_argument(u)
         case (arg_custom)
            if (present(custom_argument)) then
               if (associated(custom_argument)) then
                  argument = custom_argument(u)
               else
                  argument = u
               end if
            else
               argument = u
            end if
         case default
            argument = u
         end select
         if (trig_kind == trig_sine) then
            basis = sin(argument)
         else
            basis = cos(argument)
         end if
         values(i) = coefficients(1)
         if (size(coefficients) > 1) values(i) = values(i)+sum(coefficients(2:))*basis
      end do
   end function trigonometric_values

   function evaluate_tv_function(spec, n, status, message) result(values)
      type(tv_function_spec), intent(in) :: spec
      integer, intent(in) :: n
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      real(dp), allocatable :: values(:)

      if (present(status)) status = 0
      if (present(message)) message = 'ok'
      if (.not. allocated(spec%coefficients) .or. size(spec%coefficients) == 0) then
         allocate(values(max(0,n)))
         values = 0.0_dp
         if (present(status)) status = 1
         if (present(message)) message = 'coefficient vector is empty'
         return
      end if

      select case (spec%kind)
      case (tv_polynomial)
         values = polynomial_values(n,spec%coefficients)
      case (tv_nonlinear)
         if (.not. allocated(spec%exponents) .or. size(spec%exponents) < size(spec%coefficients)) then
            allocate(values(max(0,n)))
            values = 0.0_dp
            if (present(status)) status = 2
            if (present(message)) message = 'nonlinear exponents are missing or too short'
            return
         end if
         values = nonlinear_values(n,spec%coefficients,spec%exponents)
      case (tv_trigonometric)
         if (spec%argument_kind == arg_custom .and. associated(spec%custom_argument)) then
            values = trigonometric_values(n,spec%coefficients,spec%trig_kind,spec%argument_kind,spec%custom_argument)
         else
            values = trigonometric_values(n,spec%coefficients,spec%trig_kind,spec%argument_kind)
         end if
      case default
         allocate(values(max(0,n)))
         values = 0.0_dp
         if (present(status)) status = 3
         if (present(message)) message = 'invalid time-varying function type'
         return
      end select

      if (any(.not. ieee_is_finite(values))) then
         if (present(status)) status = 4
         if (present(message)) message = 'time-varying function produced a non-finite value'
      end if
   end function evaluate_tv_function

   pure function identity_argument(u) result(value)
      real(dp), intent(in) :: u
      real(dp) :: value
      value = u
   end function identity_argument

   pure function three_one_minus_log_argument(u) result(value)
      real(dp), intent(in) :: u
      real(dp) :: value
      value = 3.0_dp*(1.0_dp-log(max(u,tiny(1.0_dp))))
   end function three_one_minus_log_argument

end module tvgarchkf_functions
