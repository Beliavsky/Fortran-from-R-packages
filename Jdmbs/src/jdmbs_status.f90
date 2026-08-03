! SPDX-License-Identifier: GPL-2.0-or-later
module jdmbs_status
   implicit none
   integer, parameter :: jdmbs_success = 0
   integer, parameter :: jdmbs_invalid_argument = 1
   integer, parameter :: jdmbs_nonfinite_input = 2
   integer, parameter :: jdmbs_numerical_warning = 3
   public
end module jdmbs_status
