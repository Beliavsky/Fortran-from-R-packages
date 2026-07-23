! Part of the experimental modern Fortran translation of fGarch 4052.93.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original fGarch authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-or-later

program test_distributions
   use fgarch, only : dp, dnorm_fg, dstd, dged, pstd, qstd, pged, qged, &
      psnorm, qsnorm, psstd, qsstd, psged, qsged, absolute_moment, &
      dist_norm, dist_std, dist_ged
   implicit none

   real(dp), parameter :: tol = 2.0e-9_dp
   real(dp), parameter :: pvals(5) = [0.01_dp,0.10_dp,0.50_dp,0.90_dp,0.99_dp]
   real(dp) :: q(3), pback(3)
   integer :: i

   call assert_close(dnorm_fg(0.0_dp,0.0_dp,1.0_dp),0.3989422804014327_dp,tol,'normal density')
   call assert_close(dged(0.0_dp,0.0_dp,1.0_dp,2.0_dp),0.3989422804014327_dp,tol,'GED nu=2 density')
   call assert_close(dstd(0.0_dp,0.0_dp,1.0_dp,5.0_dp),0.4900701292638152_dp,2.0e-9_dp,'std density')

   do i = 1, size(pvals)
      call assert_close(pstd(qstd(pvals(i),0.0_dp,1.0_dp,7.0_dp),0.0_dp,1.0_dp,7.0_dp), &
                        pvals(i),2.0e-10_dp,'std cdf/quantile')
      call assert_close(pged(qged(pvals(i),0.0_dp,1.0_dp,1.4_dp),0.0_dp,1.0_dp,1.4_dp), &
                        pvals(i),2.0e-10_dp,'GED cdf/quantile')
   end do

   q = [qsnorm(0.49_dp,0.0_dp,1.0_dp,1.5_dp), &
        qsnorm(0.50_dp,0.0_dp,1.0_dp,1.5_dp), &
        qsnorm(0.51_dp,0.0_dp,1.0_dp,1.5_dp)]
   if (any(q(2:3) <= q(1:2))) error stop 'skew-normal quantile is not monotone'
   pback = psnorm(q,0.0_dp,1.0_dp,1.5_dp)
   if (maxval(abs(pback-[0.49_dp,0.50_dp,0.51_dp])) > 5.0e-9_dp) error stop 'skew-normal round trip'

   q = [qsstd(0.49_dp,0.0_dp,1.0_dp,6.0_dp,1.5_dp), &
        qsstd(0.50_dp,0.0_dp,1.0_dp,6.0_dp,1.5_dp), &
        qsstd(0.51_dp,0.0_dp,1.0_dp,6.0_dp,1.5_dp)]
   if (any(q(2:3) <= q(1:2))) error stop 'skew-t quantile is not monotone'
   pback = psstd(q,0.0_dp,1.0_dp,6.0_dp,1.5_dp)
   if (maxval(abs(pback-[0.49_dp,0.50_dp,0.51_dp])) > 5.0e-9_dp) error stop 'skew-t round trip'

   q = [qsged(0.49_dp,0.0_dp,1.0_dp,1.5_dp,0.8_dp), &
        qsged(0.50_dp,0.0_dp,1.0_dp,1.5_dp,0.8_dp), &
        qsged(0.51_dp,0.0_dp,1.0_dp,1.5_dp,0.8_dp)]
   if (any(q(2:3) <= q(1:2))) error stop 'skew-GED quantile is not monotone'
   pback = psged(q,0.0_dp,1.0_dp,1.5_dp,0.8_dp)
   if (maxval(abs(pback-[0.49_dp,0.50_dp,0.51_dp])) > 5.0e-9_dp) error stop 'skew-GED round trip'

   call assert_close(absolute_moment(2.0_dp,dist_norm,2.0_dp),1.0_dp,1.0e-12_dp,'normal variance')
   call assert_close(absolute_moment(2.0_dp,dist_std,7.0_dp),1.0_dp,1.0e-12_dp,'std variance')
   call assert_close(absolute_moment(2.0_dp,dist_ged,1.5_dp),1.0_dp,1.0e-12_dp,'GED variance')

   print '(a)', 'distribution tests passed'

contains

   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      if (abs(actual-expected) > tolerance) then
         print '(a,2es24.14)', trim(label)//': ',actual,expected
         error stop 'assert_close failed'
      end if
   end subroutine assert_close

end program test_distributions
