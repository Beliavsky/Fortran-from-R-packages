! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from highfrequency 1.0.2 by Kris Boudt, Jonathan Cornelissen,
! Scott Payseur, Onno Kleen, Emil Sjoerup, and contributors.
program test_realized
  use highfrequency, only: dp, rrvar, rskew, rkurt, rbpvar, rquar
  use highfrequency, only: rminrvar, rmedrvar, rminrquar, rmedrquar
  use highfrequency, only: preaveraged_covariance, average_realized_covariance
  use highfrequency, only: rtpquar, rqpvar, rcov, rbpcov, rsvar
  implicit none
  real(dp) :: r(6), m(6,2), cov(2,2), bp(2,2), semi(2), robust(2,2)

  r=[0.01_dp,-0.02_dp,0.015_dp,0.005_dp,-0.01_dp,0.03_dp]
  call assert_close(rrvar(r),0.00175_dp,1.0e-14_dp)
  call assert_close(rskew(r),0.7528371991317254_dp,1.0e-13_dp)
  call assert_close(rkurt(r),2.04_dp,1.0e-13_dp)
  call assert_close(rbpvar(r),0.0017435839227423352_dp,1.0e-13_dp)
  call assert_close(rquar(r),2.0825e-6_dp,1.0e-16_dp)
  call assert_close(rminrvar(r),0.0015686048845139423_dp,1.0e-13_dp)
  call assert_close(rmedrvar(r),0.0013838743444718798_dp,1.0e-13_dp)
  call assert_close(rminrquar(r),1.1410719726144559e-6_dp,1.0e-16_dp)
  call assert_close(rmedrquar(r),1.0075528397411956e-6_dp,1.0e-16_dp)
  call assert_close(rtpquar(r),1.3247041162738416e-6_dp,1.0e-16_dp)
  call assert_close(rqpvar(r),1.554462693171574e-6_dp,1.0e-16_dp)

  m(:,1)=r
  m(:,2)=2.0_dp*r
  cov=rcov(m)
  call assert_close(cov(1,1),0.00175_dp,1.0e-14_dp)
  call assert_close(cov(1,2),0.0035_dp,1.0e-14_dp)
  call assert_close(cov(2,2),0.007_dp,1.0e-14_dp)
  bp=rbpcov(m)
  call assert_close(bp(1,2),2.0_dp*bp(1,1),1.0e-13_dp)
  robust=preaveraged_covariance(m,3,force_psd=.true.)
  if(robust(1,1) < -1.0e-15_dp .or. robust(2,2) < -1.0e-15_dp)error stop 1
  robust=average_realized_covariance(m,2)
  if(abs(robust(1,2)-2.0_dp*robust(1,1))>1.0e-13_dp)error stop 1
  semi=rsvar(r)
  call assert_close(sum(semi),rrvar(r),1.0e-14_dp)
  print '(a)', 'test_realized: PASS'
contains
  subroutine assert_close(actual,expected,tolerance)
    real(dp),intent(in)::actual,expected,tolerance
    if(abs(actual-expected)>tolerance)then
      print '(a,3es24.16)', 'mismatch: ',actual,expected,abs(actual-expected)
      error stop 1
    end if
  end subroutine assert_close
end program test_realized
