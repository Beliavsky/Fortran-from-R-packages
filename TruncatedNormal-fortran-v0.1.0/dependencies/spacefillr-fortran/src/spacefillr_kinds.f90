module spacefillr_kinds
use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
implicit none
private
public :: int32, int64, real32, real64, int128
integer, parameter :: int128 = selected_int_kind(38)
end module spacefillr_kinds
