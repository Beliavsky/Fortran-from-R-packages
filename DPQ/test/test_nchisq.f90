program test_nchisq
   use r_compat, only: dp
   use dpq_nchisq
   use dpq_wiener
   implicit none
   real(dp)::p,q,d,pw
   integer::ifault
   call check_close(pnchisq(7.25_dp,4.5_dp,3.2_dp),0.5368520791283207_dp,3.0e-12_dp,'ncx2 cdf')
   call check_close(dnchisq_r(7.25_dp,4.5_dp,3.2_dp),0.08639923920902227_dp,3.0e-12_dp,'ncx2 pdf')
   q=qnchisq(0.73_dp,4.5_dp,3.2_dp)
   call check_close(q,9.884973807921401_dp,3.0e-11_dp,'ncx2 quantile')
   p=pnchisq(q,4.5_dp,3.2_dp)
   call check_close(p,0.73_dp,3.0e-11_dp,'ncx2 pq inversion')
   d=dnchisq_bessel(7.25_dp,4.5_dp,3.2_dp)
   call check_close(d,dnchisq_r(7.25_dp,4.5_dp,3.2_dp),2.0e-10_dp,'Bessel density')
   call check_close(pnchi1sq(3.2_dp,1.7_dp),pnchisq(3.2_dp,1.0_dp,1.7_dp),2.0e-12_dp,'df=1 special')
   pw=pchisq_w(7.25_dp,4.5_dp,3.2_dp,ifault=ifault)
   if (ifault /= 0 .or. pw < 0.0_dp .or. pw > 1.0_dp) error stop 'Wiener approximation'
   if (.not.(pnchisq_patnaik(7.25_dp,4.5_dp,3.2_dp) > 0.0_dp)) error stop 'Patnaik'
   print '(a)', 'test_nchisq: PASS'
contains
   subroutine check_close(got,want,tol,label)
      real(dp),intent(in)::got,want,tol
      character(len=*),intent(in)::label
      if (.not.(abs(got-want) <= tol*max(1.0_dp,abs(want)))) then
         write(*,'(a,2es24.15)') trim(label)//' failed: ',got,want
         error stop 1
      end if
   end subroutine
end program
