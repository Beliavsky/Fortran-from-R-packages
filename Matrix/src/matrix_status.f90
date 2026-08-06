! SPDX-License-Identifier: GPL-3.0-only
module matrix_status
   implicit none
   private
   integer, parameter, public :: matrix_success = 0
   integer, parameter, public :: matrix_err_shape = 1
   integer, parameter, public :: matrix_err_singular = 2
   integer, parameter, public :: matrix_err_not_posdef = 3
   integer, parameter, public :: matrix_err_invalid = 4
   integer, parameter, public :: matrix_err_convergence = 5
   integer, parameter, public :: matrix_err_io = 6
end module matrix_status
