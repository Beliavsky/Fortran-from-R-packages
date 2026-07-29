! SPDX-License-Identifier: GPL-3.0-or-later
module acdm_kinds
  implicit none
  private

  integer, parameter, public :: dp = kind(1.0d0)
  real(dp), parameter, public :: pi = acos(-1.0_dp)
  real(dp), parameter, public :: sqrt_two_pi = sqrt(2.0_dp * pi)
  real(dp), parameter, public :: tiny_pos = 1.0e-12_dp
  real(dp), parameter, public :: huge_penalty = 1.0e100_dp

  integer, parameter, public :: ACDM_SUCCESS = 0
  integer, parameter, public :: ACDM_BAD_INPUT = 1
  integer, parameter, public :: ACDM_BAD_PARAMETER = 2
  integer, parameter, public :: ACDM_NUMERIC_FAILURE = 3
  integer, parameter, public :: ACDM_NOT_CONVERGED = 4
  integer, parameter, public :: ACDM_SINGULAR = 5

end module acdm_kinds
