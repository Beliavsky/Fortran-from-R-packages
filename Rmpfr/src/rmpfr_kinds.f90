module rmpfr_kinds
  use iso_fortran_env, only: real64, int64
  implicit none
  private

  integer, parameter, public :: dp = real64
  integer, parameter, public :: i64 = int64
end module rmpfr_kinds
