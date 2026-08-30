! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of gbm3 3.0.3; see NOTICE.md and upstream/.
! gbm3-fortran - modern Fortran translation of gbm3 computational code.
! Upstream gbm3 is GPL-2.0-or-later. See LICENSE and NOTICE.md.
module gbm3_kinds
   use iso_fortran_env, only : real64
   implicit none
   private
   integer, parameter, public :: dp = real64
end module gbm3_kinds
