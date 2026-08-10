! Matrix-fortran: modern Fortran translation of computational ideas from Matrix.
! Copyright (C) 2026 Fortran translation contributors.
! SPDX-License-Identifier: GPL-3.0-only
module matrix_kinds
   implicit none
   private
   integer, parameter, public :: dp = kind(1.0d0)
   real(dp), parameter, public :: matrix_eps = epsilon(1.0_dp)
end module matrix_kinds
