! SPDX-License-Identifier: MIT
! Modern Fortran translation of R4GoodPersonalFinances computational code.
module r4gpf_status
  implicit none
  private
  integer, parameter, public :: r4gpf_success = 0
  integer, parameter, public :: r4gpf_invalid_argument = 1
  integer, parameter, public :: r4gpf_dimension_error = 2
  integer, parameter, public :: r4gpf_numerical_error = 3
  integer, parameter, public :: r4gpf_not_converged = 4
  integer, parameter, public :: r4gpf_allocation_error = 5
end module r4gpf_status
