! Part of the experimental modern Fortran translation of tseries 0.10-62.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original tseries authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only

module tseries_kinds
   implicit none
   private
   public :: dp
   integer, parameter :: dp = kind(1.0d0)
end module tseries_kinds
