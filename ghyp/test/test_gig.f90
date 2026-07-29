! SPDX-License-Identifier: GPL-2.0-or-later
program test_gig
   use ghyp
   implicit none
   real(dp) :: x, p, q, m, v

   call assert_close(bessel_k(0.7_dp,1.3_dp),0.3214020154044203_dp,2.0e-12_dp,'bessel K')
   x = dgig(1.2_dp,0.7_dp,1.4_dp,2.3_dp)
   call assert_close(x,0.4816458319533276_dp,2.0e-11_dp,'GIG density')
   p = pgig(1.2_dp,0.7_dp,1.4_dp,2.3_dp)
   call assert_close(p,0.5564127877672688_dp,2.0e-9_dp,'GIG CDF')
   q = qgig(p,0.7_dp,1.4_dp,2.3_dp)
   call assert_close(q,1.2_dp,2.0e-8_dp,'GIG quantile')
   m = gig_mean(0.7_dp,1.4_dp,2.3_dp)
   v = gig_variance(0.7_dp,1.4_dp,2.3_dp)
   call assert_close(m,1.3214917931119539_dp,2.0e-11_dp,'GIG mean')
   call assert_close(v,0.8158647001206434_dp,4.0e-10_dp,'GIG variance')
   call assert_close(gig_mean(-3.0_dp,4.0_dp,0.0_dp),1.0_dp,1.0e-13_dp,'inverse gamma mean')
   call assert_close(gig_mean(2.0_dp,0.0_dp,4.0_dp),1.0_dp,1.0e-13_dp,'gamma mean')
   print '(a)', 'test_gig: PASS'
contains
   subroutine assert_close(actual,expected,tol,label)
      real(dp),intent(in)::actual,expected,tol
      character(len=*),intent(in)::label
      if(abs(actual-expected)>tol*(1.0_dp+abs(expected)))then
         write(*,'(a,3es24.16)')trim(label)//' mismatch: ',actual,expected,abs(actual-expected)
         error stop 1
      end if
   end subroutine assert_close
end program test_gig
