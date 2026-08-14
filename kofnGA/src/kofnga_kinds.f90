module kofnga_kinds
  use, intrinsic :: iso_fortran_env, only : int32, int64, real64
  implicit none
  private
  public :: i32, i64, dp
  integer, parameter :: i32 = int32
  integer, parameter :: i64 = int64
  integer, parameter :: dp = real64
end module kofnga_kinds
