! Numerical kinds for the modern Fortran translation of jomo.
! Upstream jomo 2.7-6 by Matteo Quartagno and James Carpenter; License: GPL-2.
! Modern Fortran translation, 2026. Distributed under GPL-2.0-only.
module jomo_kinds
   use iso_fortran_env, only : real64, int64
   implicit none
   private
   integer, parameter, public :: dp = real64
   integer, parameter, public :: i8 = int64
end module jomo_kinds
