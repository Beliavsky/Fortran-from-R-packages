! SPDX-License-Identifier: GPL-3.0-only
program test_core_statistics
  use mixedindtests
  implicit none
  real(dp) :: x(5,3)
  type(prepared_data_result) :: pd
  type(pair_dependence_result) :: pair
  type(sn_result) :: out
  type(moebius_result) :: mo

  x = reshape([1._dp,2._dp,1._dp,3._dp,2._dp, &
               3._dp,1._dp,2._dp,2._dp,1._dp, &
               0._dp,1._dp,0._dp,1._dp,2._dp],[5,3])

  pd = preparedata(x(:,1))
  call check(pd%status == mixedind_success, 'prepare status')
  call check(size(pd%values) == 3, 'unique count')
  call check(maxval(abs(pd%values-[1._dp,2._dp,3._dp])) < 1.0e-14_dp, 'unique values')
  call check(maxval(abs(pd%pdf-[0.4_dp,0.4_dp,0.2_dp])) < 1.0e-14_dp, 'empirical pdf')

  pair = stat_dep(x(:,1),x(:,2))
  call close(pair%tau,-0.24_dp,1.0e-14_dp,'tau')
  call close(pair%rho,-0.4722222222222222_dp,1.0e-14_dp,'rho')
  call close(pair%scale,0.288_dp,1.0e-14_dp,'pair scale')

  out = Sn_A(x,3)
  call close(out%sn,0.013224296296296335_dp,1.0e-14_dp,'nonserial Sn')
  call check(all(out%cardinality == [2,2,2,3]), 'nonserial cardinalities')
  call check(maxval(abs(out%stats-[0.016177777777777774_dp, &
    0.030044444444444447_dp,0.026133333333333342_dp, &
    0.0006542222222222222_dp])) < 1.0e-14_dp, 'nonserial statistics')
  call close(out%sn_multiplier(1,1),0.25502814814814823_dp,1.0e-14_dp,'J(1,1)')
  call close(out%sn_multiplier(2,3),0.14483555555555555_dp,1.0e-14_dp,'J(2,3)')
  call close(out%multiplier(1,1,1),0.0256_dp,1.0e-14_dp,'M(1,1,1)')

  mo = EstDepMoebius(x,3)
  call check(mo%status == mixedind_success, 'Moebius status')
  call close(mo%spearman(1),-0.4722222222222223_dp,1.0e-13_dp,'Moebius Spearman')
  call close(mo%vdw(2),0.6758919955595887_dp,1.0e-13_dp,'Moebius VDW')
  call close(mo%savage(4),-0.3384655344496637_dp,1.0e-13_dp,'Moebius Savage')

  print '(a)', 'test_core_statistics: PASS'
contains
  subroutine check(ok,name)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: name
    if (.not. ok) then
      print '(a)', 'FAIL: '//trim(name)
      error stop 1
    end if
  end subroutine check
  subroutine close(x0,y0,tol,name)
    real(dp), intent(in) :: x0,y0,tol
    character(len=*), intent(in) :: name
    call check(abs(x0-y0) <= tol,name)
  end subroutine close
end program test_core_statistics
