! SPDX-License-Identifier: LGPL-3.0-or-later
program test_distributions
   use sharper, only: dp, student_t_cdf, f_cdf, nct_cdf, nct_pdf
   use sharper, only: ncf_cdf, ncf_pdf, nct_quantile, ncf_quantile
   use sharper, only: psr, qsr, psropt, qsropt, plambdap, qlambdap
   implicit none
   real(dp) :: p, q

   call assert_close(student_t_cdf(1.2_dp,10.0_dp),0.8711018496378472_dp,2.0e-12_dp)
   call assert_close(f_cdf(2.0_dp,4.0_dp,20.0_dp),0.8666523497275129_dp,2.0e-12_dp)
   call assert_close(nct_cdf(1.2_dp,10.0_dp,0.7_dp),0.6751651489644870_dp,3.0e-11_dp)
   call assert_close(nct_pdf(1.2_dp,10.0_dp,0.7_dp),0.3300988958641937_dp,3.0e-11_dp)
   call assert_close(ncf_cdf(2.0_dp,4.0_dp,20.0_dp,3.0_dp),0.6254295719853862_dp,3.0e-12_dp)
   call assert_close(ncf_pdf(2.0_dp,4.0_dp,20.0_dp,3.0_dp),0.2646774318185070_dp,3.0e-12_dp)

   p = 0.83_dp
   q = nct_quantile(p,12.0_dp,0.4_dp)
   call assert_close(nct_cdf(q,12.0_dp,0.4_dp),p,2.0e-10_dp)
   q = ncf_quantile(p,5.0_dp,30.0_dp,2.5_dp)
   call assert_close(ncf_cdf(q,5.0_dp,30.0_dp,2.5_dp),p,2.0e-10_dp)

   q = qsr(0.73_dp,40.0_dp,0.4_dp,252.0_dp)
   call assert_close(psr(q,40.0_dp,0.4_dp,252.0_dp),0.73_dp,3.0e-9_dp)
   q = qsropt(0.62_dp,4.0_dp,80.0_dp,0.8_dp,12.0_dp)
   call assert_close(psropt(q,4.0_dp,80.0_dp,0.8_dp,12.0_dp),0.62_dp,3.0e-10_dp)

   q = qlambdap(0.8_dp,20.0_dp,1.4_dp)
   call assert_close(plambdap(q,20.0_dp,1.4_dp),0.8_dp,3.0e-9_dp)

   print '(a)', 'test_distributions: PASS'
contains
   subroutine assert_close(actual, expected, tolerance)
      real(dp), intent(in) :: actual, expected, tolerance
      if (abs(actual-expected) > tolerance) then
         print '(a,3es24.16)', 'assert_close failed: ',actual,expected,tolerance
         error stop 1
      end if
   end subroutine assert_close
end program test_distributions
