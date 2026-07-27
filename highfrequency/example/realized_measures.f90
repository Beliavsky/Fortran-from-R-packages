! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from highfrequency 1.0.2 by Kris Boudt, Jonathan Cornelissen,
! Scott Payseur, Onno Kleen, Emil Sjoerup, and contributors.
program realized_measures
  use highfrequency, only: dp, rcov, rbpcov, rkernelcov, rthresholdcov
  implicit none
  real(dp) :: returns(8,2), cov(2,2)
  integer :: i
  do i=1,8
    returns(i,1)=0.01_dp*sin(0.7_dp*real(i,dp))
    returns(i,2)=0.008_dp*sin(0.7_dp*real(i,dp)+0.2_dp)
  end do
  cov=rcov(returns)
  print '(a,4es14.5)', 'realized covariance: ',cov
  cov=rbpcov(returns)
  print '(a,4es14.5)', 'bipower covariance: ',cov
  cov=rkernelcov(returns,2,'parzen',force_psd=.true.)
  print '(a,4es14.5)', 'kernel covariance: ',cov
  cov=rthresholdcov(returns)
  print '(a,4es14.5)', 'threshold covariance: ',cov
end program realized_measures
