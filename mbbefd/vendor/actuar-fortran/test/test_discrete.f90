! SPDX-License-Identifier: GPL-2.0-or-later
program test_discrete
  use actuar, only : dp, dztpois, pztpois, dztgeom, dztbinom, &
    dztnbinom, dlogarithmic, dzmpois, pzmpois, dpoisinvgauss
  implicit none

  call assert_close(dztpois(3.0_dp,2.5_dp),0.2328788797756563_dp,1.0e-13_dp)
  call assert_close(pztpois(3.0_dp,2.5_dp),0.7358972600910740_dp,1.0e-13_dp)
  call assert_close(dztgeom(4.0_dp,0.35_dp),0.09611875_dp,1.0e-14_dp)
  call assert_close(dztbinom(2.0_dp,6,0.3_dp),0.3673538081783779_dp,1.0e-13_dp)
  call assert_close(dztnbinom(3.0_dp,2.5_dp,0.4_dp), &
    0.1595903195372997_dp,2.0e-13_dp)
  call assert_close(dlogarithmic(3.0_dp,0.6_dp), &
    0.0785776800914850_dp,1.0e-14_dp)
  call assert_close(dzmpois(0.0_dp,2.0_dp,0.25_dp),0.25_dp,1.0e-15_dp)
  call assert_close(pzmpois(0.0_dp,2.0_dp,0.25_dp),0.25_dp,1.0e-15_dp)
  call assert_close(sum_pig(4,3.0_dp,0.2_dp), &
    sum_pig_reference(),2.0e-13_dp)

  print '(a)', 'test_discrete: PASS'
contains
  function sum_pig(n,mu,phi) result(s)
    integer, intent(in) :: n
    real(dp), intent(in) :: mu,phi
    real(dp) :: s
    integer :: k
    s=0.0_dp
    do k=0,n
      s=s+dpoisinvgauss(real(k,dp),mu,phi)
    end do
  end function sum_pig

  pure function sum_pig_reference() result(s)
    real(dp) :: s,p0,p1,p2,p3,p4,a,b,mu,phi
    mu=3.0_dp; phi=0.2_dp
    p0=exp((1.0_dp-sqrt(1.0_dp+2.0_dp*phi*mu*mu))/(phi*mu))
    p1=mu*p0/sqrt(1.0_dp+2.0_dp*phi*mu*mu)
    a=1.0_dp/(1.0_dp+1.0_dp/(2.0_dp*phi*mu*mu))
    b=mu*mu/(1.0_dp+2.0_dp*phi*mu*mu)
    p2=a*(1.0_dp-1.5_dp/2.0_dp)*p1+b*p0/2.0_dp
    p3=a*(1.0_dp-1.5_dp/3.0_dp)*p2+b*p1/6.0_dp
    p4=a*(1.0_dp-1.5_dp/4.0_dp)*p3+b*p2/12.0_dp
    s=p0+p1+p2+p3+p4
  end function sum_pig_reference

  subroutine assert_close(actual,expected,tol)
    real(dp), intent(in) :: actual,expected,tol
    if(abs(actual-expected)>tol*max(1.0_dp,abs(expected))) then
      print '(a,3es24.15)', 'mismatch: ',actual,expected,abs(actual-expected)
      error stop 1
    end if
  end subroutine assert_close
end program test_discrete
