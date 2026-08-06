! Modern Fortran translation of the computational core of tvGarchKF.
! Copyright (C) 2026 OpenAI
! SPDX-License-Identifier: GPL-3.0-or-later
module tvgarchkf_types
   use fgarch_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: tv_polynomial = 1
   integer, parameter, public :: tv_nonlinear = 2
   integer, parameter, public :: tv_trigonometric = 3
   integer, parameter, public :: trig_sine = 1
   integer, parameter, public :: trig_cosine = 2
   integer, parameter, public :: arg_identity = 1
   integer, parameter, public :: arg_three_one_minus_log = 2
   integer, parameter, public :: arg_custom = 3

   abstract interface
      pure function tv_argument_function(u) result(value)
         import :: dp
         real(dp), intent(in) :: u
         real(dp) :: value
      end function tv_argument_function
   end interface

   type, public :: tv_function_spec
      integer :: kind = tv_polynomial
      integer :: trig_kind = trig_cosine
      integer :: argument_kind = arg_identity
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: exponents(:)
      procedure(tv_argument_function), pointer, nopass :: custom_argument => null()
   end type tv_function_spec

   type, public :: tvgarch_spec
      type(tv_function_spec) :: omega
      type(tv_function_spec) :: alpha
      type(tv_function_spec) :: beta
   end type tvgarch_spec

   type, public :: tvgarch_filter_result
      integer :: status = 1
      character(len=160) :: message = 'not evaluated'
      real(dp) :: criterion = huge(1.0_dp)
      real(dp), allocatable :: omega(:)
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: state(:)
      real(dp), allocatable :: state_variance(:)
      real(dp), allocatable :: mse(:)
      real(dp), allocatable :: gain(:)
      real(dp), allocatable :: conditional_variance(:)
      real(dp), allocatable :: sigma(:)
   end type tvgarch_filter_result

   type, public :: tvgarch_fit_result
      integer :: status = 1
      character(len=160) :: message = 'not fitted'
      integer :: iterations = 0
      integer :: evaluations = 0
      real(dp) :: criterion = huge(1.0_dp)
      real(dp), allocatable :: parameters(:)
      real(dp), allocatable :: rounded_parameters(:)
      type(tvgarch_spec) :: spec
      type(tvgarch_filter_result) :: filter
   end type tvgarch_fit_result

   type, public :: tvgarch_simulation_result
      integer :: status = 1
      character(len=160) :: message = 'not simulated'
      real(dp), allocatable :: returns(:)
      real(dp), allocatable :: variance(:)
      real(dp), allocatable :: omega(:)
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: beta(:)
   end type tvgarch_simulation_result

   type, public :: tv_parameter_result
      integer :: status = 1
      character(len=160) :: message = 'not estimated'
      real(dp), allocatable :: midpoint(:)
      real(dp), allocatable :: omega(:)
      real(dp), allocatable :: alpha(:)
      real(dp), allocatable :: beta(:)
      real(dp) :: global_omega = 0.0_dp
      real(dp) :: global_alpha = 0.0_dp
      real(dp) :: global_beta = 0.0_dp
   end type tv_parameter_result

   public :: make_tv_function, make_tvgarch_spec, tv_argument_function

contains

   function make_tv_function(coefficients, kind, exponents, trig_kind, argument_kind) result(spec)
      real(dp), intent(in) :: coefficients(:)
      integer, intent(in), optional :: kind, trig_kind, argument_kind
      real(dp), intent(in), optional :: exponents(:)
      type(tv_function_spec) :: spec

      allocate(spec%coefficients(size(coefficients)))
      spec%coefficients = coefficients
      if (present(kind)) spec%kind = kind
      if (present(trig_kind)) spec%trig_kind = trig_kind
      if (present(argument_kind)) spec%argument_kind = argument_kind
      if (present(exponents)) then
         allocate(spec%exponents(size(exponents)))
         spec%exponents = exponents
      end if
   end function make_tv_function

   function make_tvgarch_spec(omega, alpha, beta) result(spec)
      type(tv_function_spec), intent(in) :: omega, alpha, beta
      type(tvgarch_spec) :: spec

      spec%omega = omega
      spec%alpha = alpha
      spec%beta = beta
   end function make_tvgarch_spec

end module tvgarchkf_types
