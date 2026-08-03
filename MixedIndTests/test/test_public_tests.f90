! SPDX-License-Identifier: GPL-3.0-only
program test_public_tests
  use iso_fortran_env, only : int64
  use mixedindtests
  implicit none
  real(dp) :: x(12,3), s(12)
  integer :: i
  type(dependence_result) :: dep
  type(serial_dependence_result) :: sdep
  type(copula_test_result) :: ind,ser,multi
  type(moebius_result) :: mo

  do i=1,12
    x(i,1)=real(mod(3*i+1,7),dp)
    x(i,2)=real(mod(5*i+2,11),dp)
    x(i,3)=sin(real(i,dp))
  end do
  s=x(:,1)

  dep=EstDep(x)
  call check(dep%status==mixedind_success,'EstDep status')
  call check(all(abs(diagonal(dep%tau)-1.0_dp)<1.0e-14_dp),'tau diagonal')
  call check(dep%p_lb_tau>=0.0_dp .and. dep%p_lb_tau<=100.0_dp,'LB p-value')

  sdep=EstDepSerial(s,3)
  call check(sdep%status==mixedind_success,'EstDepSerial status')
  call check(size(sdep%tau)==3,'serial lag count')

  ind=TestIndCopula(x,3,24,101_int64)
  call check(ind%status==mixedind_success,'independence test status')
  call probabilities(ind)
  ser=TestIndSerCopula(s,4,3,24,102_int64)
  call check(ser%status==mixedind_success,'serial test status')
  call probabilities(ser)
  multi=TestIndSerCopulaMulti(x(:,1:2),3,2,16,103_int64)
  call check(multi%status==mixedind_success,'multivariate serial status')
  call probabilities(multi)

  mo=EstDepSerialMoebius(s,4,3)
  call check(mo%status==mixedind_success,'serial Moebius status')
  call check(size(mo%spearman)==6,'serial Moebius subset count')
  call check(mo%p_ln_spearman>=0.0_dp .and. mo%p_ln_spearman<=100.0_dp, &
    'serial Moebius p-value')

  print '(a)', 'test_public_tests: PASS'
contains
  function diagonal(a) result(d)
    real(dp), intent(in) :: a(:,:)
    real(dp) :: d(min(size(a,1),size(a,2)))
    integer :: j
    do j=1,size(d); d(j)=a(j,j); end do
  end function diagonal
  subroutine probabilities(out)
    type(copula_test_result), intent(in) :: out
    call check(all(out%p_cvm>=0.0_dp .and. out%p_cvm<=100.0_dp),'CVM p-values')
    call check(out%p_sn>=0.0_dp .and. out%p_sn<=100.0_dp,'Sn p-value')
    call check(out%p_tn>=0.0_dp .and. out%p_tn<=100.0_dp,'Tn p-value')
  end subroutine probabilities
  subroutine check(ok,name)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: name
    if (.not. ok) then
      print '(a)', 'FAIL: '//trim(name)
      error stop 1
    end if
  end subroutine check
end program test_public_tests
