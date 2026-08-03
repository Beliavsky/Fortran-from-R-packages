! SPDX-License-Identifier: GPL-2.0-or-later
module infoset_status
  implicit none
  private
  integer, parameter, public :: infoset_success = 0
  integer, parameter, public :: infoset_invalid_argument = 1
  integer, parameter, public :: infoset_insufficient_data = 2
  integer, parameter, public :: infoset_no_split = 3
  integer, parameter, public :: infoset_not_converged = 4
  integer, parameter, public :: infoset_numerical_error = 5
  integer, parameter, public :: infoset_qp_error = 6
end module infoset_status
