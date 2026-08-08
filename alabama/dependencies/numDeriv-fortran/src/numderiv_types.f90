! SPDX-License-Identifier: GPL-2.0-or-later
module numderiv_types
   use numderiv_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: nd_success = 0
   integer, parameter, public :: nd_invalid_argument = 1
   integer, parameter, public :: nd_nonfinite_value = 2
   integer, parameter, public :: nd_shape_mismatch = 3

   type, public :: deriv_options
      real(dp) :: eps = 1.0e-4_dp
      real(dp) :: d = 1.0e-4_dp
      real(dp) :: zero_tol = sqrt(epsilon(1.0_dp) / 7.0e-7_dp)
      integer :: r = 4
      real(dp) :: v = 2.0_dp
      logical :: show_details = .false.
   end type deriv_options

   type, public :: gend_result
      real(dp), allocatable :: dmat(:, :)
      real(dp), allocatable :: f0(:)
      real(dp), allocatable :: x(:)
      integer :: p = 0
      real(dp) :: d = 0.0_dp
      type(deriv_options) :: options
      integer :: status = nd_success
      character(:), allocatable :: message
   end type gend_result

   public :: first_deriv_options, hessian_options

contains

   pure function first_deriv_options() result(options)
      type(deriv_options) :: options
      options = deriv_options()
   end function first_deriv_options

   pure function hessian_options() result(options)
      type(deriv_options) :: options
      options = deriv_options(d=0.1_dp)
   end function hessian_options

end module numderiv_types
