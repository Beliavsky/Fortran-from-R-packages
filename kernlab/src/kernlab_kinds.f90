! SPDX-License-Identifier: GPL-2.0-only
module kernlab_kinds
  implicit none
  private
  integer, parameter, public :: dp = kind(1.0d0)
  integer, parameter, public :: KL_SUCCESS = 0
  integer, parameter, public :: KL_INVALID_ARGUMENT = 1
  integer, parameter, public :: KL_SINGULAR = 2
  integer, parameter, public :: KL_NOT_CONVERGED = 3
  integer, parameter, public :: KL_NUMERICAL_ERROR = 4
end module kernlab_kinds
