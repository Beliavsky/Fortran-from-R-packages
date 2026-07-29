! SPDX-License-Identifier: GPL-2.0-or-later
program test_fitting
   use ghyp
   implicit none
   real(dp) :: data(10), mvdata(6,2)
   type(fit_result) :: uvfit, mvfit, quickfit
   type(likelihood_ratio_result) :: lr

   data=[-1.0_dp,-0.5_dp,-0.2_dp,0.0_dp,0.1_dp,0.3_dp,0.6_dp,0.9_dp,1.2_dp,1.6_dp]
   uvfit=fit_gaussian_uv(data)
   call assert_true(uvfit%ok .and. uvfit%converged,'Gaussian univariate fit')
   call assert_close(uvfit%model%mu(1),0.3_dp,2.0e-13_dp,'fitted mean')
   call assert_true(uvfit%aic<huge(1.0_dp),'finite AIC')

   quickfit=fit_ghyp_uv(data,'nig',max_iter=8)
   call assert_true(quickfit%ok,'NIG fitting path')

   mvdata=reshape([ &
      -0.5_dp,-0.2_dp,0.0_dp,0.2_dp,0.4_dp,0.8_dp, &
       0.3_dp,0.0_dp,-0.1_dp,0.1_dp,0.5_dp,0.9_dp],[6,2])
   mvfit=fit_gaussian_mv(mvdata)
   call assert_true(mvfit%ok .and. mvfit%model%dimension()==2,'Gaussian multivariate fit')
   lr=likelihood_ratio_test(uvfit,uvfit)
   call assert_true(.not.lr%ok,'LR degrees-of-freedom validation')
   print '(a)', 'test_fitting: PASS'
contains
   subroutine assert_close(actual,expected,tol,label)
      real(dp),intent(in)::actual,expected,tol
      character(len=*),intent(in)::label
      if(abs(actual-expected)>tol*(1.0_dp+abs(expected)))then
         write(*,'(a,3es24.16)')trim(label)//' mismatch: ',actual,expected,abs(actual-expected)
         error stop 1
      end if
   end subroutine assert_close
   subroutine assert_true(condition,label)
      logical,intent(in)::condition
      character(len=*),intent(in)::label
      if(.not.condition)then;write(*,'(a)')trim(label)//' failed';error stop 1;end if
   end subroutine assert_true
end program test_fitting
