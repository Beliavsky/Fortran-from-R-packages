! SPDX-License-Identifier: GPL-2.0-only
! Copyright (C) 2014-2025 Genaro Sucarrat (original R/C package)
! Copyright (C) 2026 Modern Fortran translation contributors
! This file is distributed under the GNU General Public License version 2 only.
program test_utils
  use lgarch_kinds, only : dp
  use lgarch_utils, only : lag_vector, lag_matrix, diff_vector, diff_matrix
  implicit none
  real(dp) :: x(5),m(4,2)
  real(dp),allocatable :: y(:),a(:,:)
  x=[1.0_dp,2.0_dp,4.0_dp,7.0_dp,11.0_dp]
  y=lag_vector(x,2,pad=.true.,pad_value=-99.0_dp)
  call assert_close(y,[-99.0_dp,-99.0_dp,1.0_dp,2.0_dp,4.0_dp],1.0e-14_dp,"lag vector padded")
  y=lag_vector(x,2,pad=.false.)
  call assert_close(y,[1.0_dp,2.0_dp,4.0_dp],1.0e-14_dp,"lag vector unpadded")
  y=diff_vector(x,2,pad=.true.,pad_value=0.0_dp)
  call assert_close(y,[0.0_dp,0.0_dp,3.0_dp,5.0_dp,7.0_dp],1.0e-14_dp,"diff vector")
  m=reshape([1.0_dp,2.0_dp,3.0_dp,4.0_dp,10.0_dp,20.0_dp,30.0_dp,40.0_dp],[4,2])
  a=lag_matrix(m,1,pad=.true.,pad_value=0.0_dp)
  call assert_close(reshape(a,[size(a)]),[0.0_dp,1.0_dp,2.0_dp,3.0_dp,0.0_dp,10.0_dp,20.0_dp,30.0_dp],1.0e-14_dp,"lag matrix")
  a=diff_matrix(m,1,pad=.false.)
  call assert_close(reshape(a,[size(a)]),[1.0_dp,1.0_dp,1.0_dp,10.0_dp,10.0_dp,10.0_dp],1.0e-14_dp,"diff matrix")
  print '(a)', 'Utility tests passed.'
contains
  subroutine assert_close(actual,expected,tol,label)
    real(dp),intent(in)::actual(:),expected(:),tol
    character(len=*),intent(in)::label
    if(size(actual)/=size(expected) .or. maxval(abs(actual-expected))>tol) then
      print *,trim(label),actual,expected
      error stop 1
    end if
  end subroutine assert_close
end program test_utils
