module suppdists_kinds
   use, intrinsic :: iso_fortran_env, only : real64, int64
   implicit none
   private
   public :: dp, i8, pi, sqrt2, sqrt2pi
   integer, parameter :: dp = real64
   integer, parameter :: i8 = int64
   real(dp), parameter :: pi = acos(-1.0_dp)
   real(dp), parameter :: sqrt2 = sqrt(2.0_dp)
   real(dp), parameter :: sqrt2pi = sqrt(2.0_dp*pi)
end module suppdists_kinds
