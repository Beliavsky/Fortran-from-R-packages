! SPDX-License-Identifier: GPL-2.0-only
!
! Modern Fortran translation of computational methods from GCPM 1.2.2.
! Original software copyright (C) 2015 Kevin Jakob and Dr. Matthias Fischer.
! Fortran translation copyright (C) 2026.
module gcpm_kinds
   implicit none
   private

   integer, parameter, public :: dp = kind(1.0d0)
   integer, parameter, public :: name_len = 128
   integer, parameter, public :: sector_name_len = 64

end module gcpm_kinds
