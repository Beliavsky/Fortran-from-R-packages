! SPDX-License-Identifier: GPL-3.0-only
program test_validation
  use indgenerrors
  implicit none
  real(dp) :: x(5) = [1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
  real(dp) :: y(4) = [1.0_dp,2.0_dp,3.0_dp,4.0_dp]
  real(dp) :: constant(5) = 1.0_dp
  type(cvm_test_result) :: a
  type(lag_test_result) :: b
  type(dependence_two_result) :: c

  a = cvm_2series(x,y,1)
  call check(a%status == indgen_invalid_argument,'unequal lengths')
  b = crosscor_2series(x,x,-1)
  call check(b%status == indgen_invalid_argument,'negative lag')
  b = crosscor_2series(constant,x,1)
  call check(b%status == indgen_numerical_error,'constant-series correlation')
  c = crossdep_2series(constant,x,1)
  call check(c%status == indgen_numerical_error,'constant-series dependence')
  print '(a)', 'test_validation: PASS'

contains

  subroutine check(condition,message)
    logical, intent(in) :: condition
    character(*), intent(in) :: message
    if (.not. condition) then
      print '(a)', 'FAIL: '//message
      error stop 1
    end if
  end subroutine check

end program test_validation
