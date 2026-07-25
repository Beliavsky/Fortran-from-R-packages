! SPDX-License-Identifier: GPL-2.0-or-later
!
! Computational translation of MSGARCH, copyright (C) MSGARCH authors.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 2 or later.
module msgarch_kinds
   implicit none
   private
   integer, parameter, public :: dp = kind(1.0d0)
   real(dp), parameter, public :: pi = acos(-1.0_dp)
   real(dp), parameter, public :: sqrt_two = sqrt(2.0_dp)
end module msgarch_kinds
