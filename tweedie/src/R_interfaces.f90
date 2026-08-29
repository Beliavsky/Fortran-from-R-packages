! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from / supporting R package tweedie 3.1.0 by Peter K. Dunn.
module R_interfaces
use, intrinsic :: iso_c_binding, only: c_int, c_double
implicit none
private
public :: intpr, dblepr
contains
subroutine intpr(label, nchar, x, n)
character(len=*), intent(in) :: label
integer, intent(in) :: nchar
integer(c_int), intent(in) :: x
integer, intent(in) :: n
if (nchar < -huge(0) .or. n < -huge(0)) continue
write(*,'(a,1x,i0)') trim(label), x
end subroutine intpr

subroutine dblepr(label, nchar, x, n)
character(len=*), intent(in) :: label
integer, intent(in) :: nchar
real(c_double), intent(in) :: x
integer, intent(in) :: n
if (nchar < -huge(0) .or. n < -huge(0)) continue
write(*,'(a,1x,es24.16)') trim(label), x
end subroutine dblepr
end module R_interfaces
