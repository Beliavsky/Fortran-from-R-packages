! SPDX-License-Identifier: GPL-3.0-only
program test_serial_statistics
  use iso_fortran_env, only : int64
  use mixedindtests
  implicit none
  real(dp) :: x(5), y(5,2), xi(5)
  type(sn_result) :: out, vec
  type(pair_dependence_result) :: dep
  type(bootstrap_result) :: bs1, bs2

  x = [1._dp,2._dp,1._dp,3._dp,2._dp]
  y = reshape([1._dp,2._dp,1._dp,3._dp,2._dp, &
               0._dp,1._dp,0._dp,1._dp,2._dp],[5,2])
  xi = [-1._dp,0.5_dp,0.25_dp,1.25_dp,-0.75_dp]

  dep = stat_dep_ser(x,1)
  call check(dep%status == mixedind_success, 'serial pair status')

  out = Sn_Aserial(x,3,3)
  call close(out%sn,0.01053629629629633_dp,1.0e-14_dp,'serial Sn')
  call check(all(out%cardinality == [2,2,3]), 'serial cardinalities')
  call check(maxval(abs(out%stats-[0.01617777777777777_dp, &
    0.005866666666666663_dp,0.0006542222222222211_dp])) < 1.0e-14_dp, &
    'serial statistics')
  call close(out%sn_multiplier(1,1),0.011527111111111111_dp,1.0e-14_dp,'serial J11')
  call close(out%sn_multiplier(2,3),0.003268740740740739_dp,1.0e-14_dp,'serial J23')

  vec = Sn_AserialVec(y,3,3)
  call close(vec%sn,0.0038707685276268915_dp,1.0e-14_dp,'vector serial Sn')
  call check(maxval(abs(vec%stats-[0.016969422222222232_dp, &
    0.002254656790123457_dp,0.001695763532729769_dp])) < 1.0e-14_dp, &
    'vector serial statistics')
  call close(vec%sn_multiplier(2,3),0.0019458436888230451_dp,1.0e-14_dp, &
    'vector serial J23')

  bs1 = bootstrap(out%multiplier,out%sn_multiplier,xi)
  bs2 = bootstrap(out%multiplier,out%sn_multiplier,xi)
  call check(maxval(abs(bs1%cvm-bs2%cvm)) < tiny(1.0_dp), 'bootstrap deterministic input')
  call close(bs1%sn,bs2%sn,0.0_dp,'bootstrap Sn repeat')
  call close(Sn_serial(x,3),out%sn,1.0e-14_dp,'Sn_serial wrapper')

  print '(a)', 'test_serial_statistics: PASS'
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
end program test_serial_statistics
