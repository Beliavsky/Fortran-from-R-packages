! SPDX-License-Identifier: GPL-3.0-only
module stochvoltmb_status
  implicit none
  private
  integer, parameter, public :: sv_ok = 0
  integer, parameter, public :: sv_invalid_argument = 1
  integer, parameter, public :: sv_no_convergence = 2
  integer, parameter, public :: sv_singular = 3
  integer, parameter, public :: sv_numerical_failure = 4
end module stochvoltmb_status
