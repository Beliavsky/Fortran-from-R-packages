! SPDX-License-Identifier: GPL-2.0-or-later
!
! Modern Fortran port of RM2006 0.1.1.
! Original R package copyright (c) 2020 Carlos Trucios.
module rm2006_kinds
   implicit none
   private

   integer, parameter, public :: dp = kind(1.0d0)
end module rm2006_kinds
