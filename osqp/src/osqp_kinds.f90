! SPDX-License-Identifier: Apache-2.0
module osqp_kinds
   use, intrinsic :: iso_c_binding, only : c_int, c_double
   implicit none
   private
   integer, parameter, public :: dp = kind(1.0d0)
   integer, parameter, public :: osqp_int = c_int
   integer, parameter, public :: osqp_real = c_double
end module osqp_kinds
