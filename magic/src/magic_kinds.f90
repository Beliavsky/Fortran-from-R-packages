module magic_kinds
   use, intrinsic :: iso_fortran_env, only : int64, real64
   implicit none
   private

   integer, parameter, public :: ik = int64
   integer, parameter, public :: dp = real64
end module magic_kinds
