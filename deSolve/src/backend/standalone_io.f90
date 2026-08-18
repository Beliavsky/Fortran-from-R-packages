! SPDX-License-Identifier: GPL-2.0-or-later
! Standalone replacements for the R printing/error callbacks used by the
! deSolve-adapted legacy solver sources.
subroutine rprintf(msg)
  use, intrinsic :: iso_fortran_env, only : error_unit
  implicit none
  character(len=*), intent(in) :: msg
  integer :: n
  n = index(msg, achar(0))
  if (n == 0) n = len_trim(msg) + 1
  if (n > 1) write(error_unit, '(a)') msg(:n-1)
end subroutine rprintf

subroutine rprintfi1(msg, i1)
  use, intrinsic :: iso_fortran_env, only : error_unit
  implicit none
  character(len=*), intent(in) :: msg
  integer, intent(in) :: i1
  write(error_unit, '(a,1x,i0)') trim_nul(msg), i1
contains
  pure function trim_nul(s) result(out)
    character(len=*), intent(in) :: s
    character(len=:), allocatable :: out
    integer :: n
    n = index(s, achar(0)); if (n == 0) n = len_trim(s) + 1
    out = s(:max(0,n-1))
  end function trim_nul
end subroutine rprintfi1

subroutine rprintfi2(msg, i1, i2)
  use, intrinsic :: iso_fortran_env, only : error_unit
  implicit none
  character(len=*), intent(in) :: msg
  integer, intent(in) :: i1, i2
  write(error_unit, '(a,1x,i0,1x,i0)') trim_nul(msg), i1, i2
contains
  pure function trim_nul(s) result(out)
    character(len=*), intent(in) :: s
    character(len=:), allocatable :: out
    integer :: n
    n = index(s, achar(0)); if (n == 0) n = len_trim(s) + 1
    out = s(:max(0,n-1))
  end function trim_nul
end subroutine rprintfi2

subroutine rprintfd1(msg, d1)
  use, intrinsic :: iso_fortran_env, only : error_unit, real64
  implicit none
  character(len=*), intent(in) :: msg
  real(real64), intent(in) :: d1
  write(error_unit, '(a,1x,es24.16)') trim_nul(msg), d1
contains
  pure function trim_nul(s) result(out)
    character(len=*), intent(in) :: s
    character(len=:), allocatable :: out
    integer :: n
    n = index(s, achar(0)); if (n == 0) n = len_trim(s) + 1
    out = s(:max(0,n-1))
  end function trim_nul
end subroutine rprintfd1

subroutine rprintfd2(msg, d1, d2)
  use, intrinsic :: iso_fortran_env, only : error_unit, real64
  implicit none
  character(len=*), intent(in) :: msg
  real(real64), intent(in) :: d1, d2
  write(error_unit, '(a,2(1x,es24.16))') trim_nul(msg), d1, d2
contains
  pure function trim_nul(s) result(out)
    character(len=*), intent(in) :: s
    character(len=:), allocatable :: out
    integer :: n
    n = index(s, achar(0)); if (n == 0) n = len_trim(s) + 1
    out = s(:max(0,n-1))
  end function trim_nul
end subroutine rprintfd2

subroutine dblepr(msg, nchar, x, n)
  use, intrinsic :: iso_fortran_env, only : error_unit, real64
  implicit none
  character(len=*), intent(in) :: msg
  integer, intent(in) :: nchar, n
  real(real64), intent(in) :: x(*)
  integer :: i
  if (nchar < -huge(1)) return
  write(error_unit, '(a)') trim(msg)
  if (n > 0) write(error_unit, '(*(1x,es24.16))') (x(i), i=1,n)
end subroutine dblepr

subroutine intpr(msg, nchar, x, n)
  use, intrinsic :: iso_fortran_env, only : error_unit
  implicit none
  character(len=*), intent(in) :: msg
  integer, intent(in) :: nchar, n, x(*)
  integer :: i
  if (nchar < -huge(1)) return
  write(error_unit, '(a)') trim(msg)
  if (n > 0) write(error_unit, '(*(1x,i0))') (x(i), i=1,n)
end subroutine intpr

subroutine rexit(msg)
  implicit none
  character(len=*), intent(in) :: msg
  error stop trim(msg)
end subroutine rexit

subroutine rprintfid(msg, i1, d1)
  use, intrinsic :: iso_fortran_env, only : error_unit, real64
  implicit none
  character(len=*), intent(in) :: msg
  integer, intent(in) :: i1
  real(real64), intent(in) :: d1
  write(error_unit,'(a,1x,i0,1x,es24.16)') trim(msg), i1, d1
end subroutine rprintfid

subroutine rprintfdi(msg, d1, i1)
  use, intrinsic :: iso_fortran_env, only : error_unit, real64
  implicit none
  character(len=*), intent(in) :: msg
  real(real64), intent(in) :: d1
  integer, intent(in) :: i1
  write(error_unit,'(a,1x,es24.16,1x,i0)') trim(msg), d1, i1
end subroutine rprintfdi

subroutine rprintfdid(msg, d1, i1, d2)
  use, intrinsic :: iso_fortran_env, only : error_unit, real64
  implicit none
  character(len=*), intent(in) :: msg
  real(real64), intent(in) :: d1,d2
  integer, intent(in) :: i1
  write(error_unit,'(a,1x,es24.16,1x,i0,1x,es24.16)') trim(msg), d1, i1, d2
end subroutine rprintfdid

subroutine rprintfi3(msg, i1, i2, i3)
  use, intrinsic :: iso_fortran_env, only : error_unit
  implicit none
  character(len=*), intent(in) :: msg
  integer, intent(in) :: i1,i2,i3
  write(error_unit,'(a,3(1x,i0))') trim(msg),i1,i2,i3
end subroutine rprintfi3
