! SPDX-License-Identifier: MIT
module invgamstochvol_status
   implicit none
   private

   integer, parameter, public :: invgam_success = 0
   integer, parameter, public :: invgam_invalid_argument = 1
   integer, parameter, public :: invgam_nonfinite_input = 2
   integer, parameter, public :: invgam_numerical_failure = 3
end module invgamstochvol_status
