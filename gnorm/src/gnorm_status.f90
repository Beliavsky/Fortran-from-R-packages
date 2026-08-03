! SPDX-License-Identifier: GPL-2.0-or-later
module gnorm_status
   implicit none
   private
   integer, parameter, public :: gnorm_success = 0
   integer, parameter, public :: gnorm_invalid_argument = 1
   integer, parameter, public :: gnorm_numerical_failure = 2
end module gnorm_status
