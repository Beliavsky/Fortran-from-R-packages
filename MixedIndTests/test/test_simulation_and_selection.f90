! SPDX-License-Identifier: GPL-3.0-only
program test_simulation_and_selection
  use iso_fortran_env, only : int64
  use mixedindtests
  implicit none
  real(dp), allocatable :: x(:),x2(:),u(:),q(:)
  integer :: p

  x=SimAR1Poisson([2.0_dp,0.3_dp],80,777_int64)
  x2=SimAR1Poisson([2.0_dp,0.3_dp],80,777_int64)
  call check(size(x)==80,'Poisson length')
  call check(maxval(abs(x-x2)) < tiny(1.0_dp),'Poisson reproducibility')
  call check(all(x>=0.0_dp),'Poisson nonnegative')

  q=Finv([0.1_dp,0.5_dp,0.9_dp],1)
  call check(maxval(abs(q-[0.0_dp,1.0_dp,1.0_dp])) < tiny(1.0_dp),'Bernoulli quantiles')
  q=Finv([0.1_dp,0.5_dp,0.9_dp],2)
  call check(all(q>=0.0_dp),'Poisson quantiles')
  q=Finv([0.2_dp,0.5_dp,0.8_dp],7)
  call check(minval(q(2:)-q(:2)) > -1.0e-14_dp,'Pareto quantiles increasing')

  u=SimCopulaSeries('gaussian',100,0.4_dp,seed=1001_int64)
  call unit_interval(u,'Gaussian copula')
  u=SimCopulaSeries('clayton',40,0.25_dp,seed=1002_int64)
  call unit_interval(u,'Clayton copula')
  u=SimCopulaSeries('frank',40,-0.2_dp,seed=1003_int64)
  call unit_interval(u,'Frank copula')
  u=SimCopulaSeries('gumbel',40,0.3_dp,seed=1004_int64)
  call unit_interval(u,'Gumbel copula')
  u=SimCopulaSeries('joe',30,0.2_dp,seed=1005_int64)
  call unit_interval(u,'Joe copula')
  u=SimCopulaSeries('plackett',40,0.2_dp,seed=1006_int64)
  call unit_interval(u,'Plackett copula')
  u=SimCopulaSeries('fgm',40,0.1_dp,3.0_dp,1007_int64)
  call unit_interval(u,'FGM copula')
  u=SimCopulaSeries('t',30,0.2_dp,5.0_dp,1008_int64)
  call unit_interval(u,'Student copula')

  p=select_p(x,2,4)
  call check(p>=2 .and. p<=4,'selected p range')

  print '(a)', 'test_simulation_and_selection: PASS'
contains
  subroutine unit_interval(v,name)
    real(dp), intent(in) :: v(:)
    character(len=*), intent(in) :: name
    call check(all(v>0.0_dp .and. v<1.0_dp),name)
  end subroutine unit_interval
  subroutine check(ok,name)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: name
    if (.not. ok) then
      print '(a)', 'FAIL: '//trim(name)
      error stop 1
    end if
  end subroutine check
end program test_simulation_and_selection
