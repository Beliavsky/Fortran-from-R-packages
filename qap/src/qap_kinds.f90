module qap_kinds
   use, intrinsic :: iso_fortran_env, only : real64, int64
   implicit none
   private
   public :: dp, i64
   integer, parameter :: dp = real64
   integer, parameter :: i64 = int64
end module qap_kinds
