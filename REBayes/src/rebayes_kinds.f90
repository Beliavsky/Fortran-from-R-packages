module rebayes_kinds
   use, intrinsic :: iso_fortran_env, only : real64, int32
   implicit none
   private
   public :: dp, i32, pi
   integer, parameter :: dp = real64
   integer, parameter :: i32 = int32
   real(dp), parameter :: pi = acos(-1.0_dp)
end module rebayes_kinds
