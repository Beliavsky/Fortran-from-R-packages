! SPDX-License-Identifier: GPL-2.0-or-later
program test_distributions
   use ghyp
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   type(ghyp_model_type) :: model, mv, gauss, stud
   type(moments_result) :: moments
   real(dp) :: sc(2,2), mu(2), gamma(2), x(2), value, q
   real(dp), allocatable :: sample(:,:)
   logical :: ok

   model = ghyp_ad(0.7_dp,1.8_dp,1.2_dp,[0.3_dp],[0.2_dp],reshape([1.0_dp],[1,1]))
   call assert_true(model%ok,'GH constructor')
   call assert_close(dghyp(0.5_dp,model),0.4337186952407726_dp,3.0e-11_dp,'GH density')
   call assert_close(pghyp(0.5_dp,model),0.5146695401849739_dp,3.0e-9_dp,'GH CDF')
   q = qghyp(0.6_dp,model)
   call assert_close(q,0.702523044295577_dp,2.0e-7_dp,'GH quantile')

   gauss = gaussian_uv(0.4_dp,1.3_dp)
   call assert_close(dghyp(0.2_dp,gauss),0.3032683817402127_dp,2.0e-13_dp,'Gaussian density')
   call assert_close(qghyp(0.975_dp,gauss),2.947953184772281_dp,3.0e-8_dp,'Gaussian quantile')

   stud = student_t_uv(6.0_dp,mu=0.1_dp,sigma=1.2_dp)
   value = dghyp(0.4_dp,stud)
   call assert_close(value,0.3699927340026566_dp,3.0e-11_dp,'Student density')

   mu=[0.1_dp,-0.2_dp]
   sc=reshape([1.0_dp,0.3_dp,0.3_dp,0.8_dp],[2,2])
   gamma=[0.2_dp,-0.1_dp]
   mv=ghyp_mv(0.7_dp,1.4_dp,2.3_dp,mu,sc,gamma)
   x=[0.4_dp,0.2_dp]
   call assert_close(dghyp(x,mv),0.1738361874059931_dp,4.0e-11_dp,'multivariate density')
   moments=ghyp_moments(mv)
   call assert_close(moments%mean(1),0.3642983586223908_dp,3.0e-11_dp,'mean 1')
   call assert_close(moments%mean(2),-0.3321491793111954_dp,3.0e-11_dp,'mean 2')
   call assert_close(moments%covariance(1,2),0.3801302449703228_dp,2.0e-9_dp,'covariance')

   call rghyp(100,mv,sample,ok,12345_i8)
   call assert_true(ok .and. size(sample,1)==100 .and. size(sample,2)==2,'simulation shape')
   call assert_true(all(ieee_is_finite(sample)),'finite simulation')
   print '(a)', 'test_distributions: PASS'
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
end program test_distributions
