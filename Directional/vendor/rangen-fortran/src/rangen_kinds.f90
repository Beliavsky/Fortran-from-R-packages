module rangen_kinds
    use, intrinsic :: iso_fortran_env, only : real64, int64
    implicit none
    private
    public :: dp, i8
    integer, parameter :: dp = real64
    integer, parameter :: i8 = int64
end module rangen_kinds
