! SPDX-License-Identifier: GPL-2.0-or-later
program test_risk_transform
   use ghyp
   implicit none
   type(ghyp_model_type) :: model, mv, transformed
   type(attribution_result) :: attrib
   type(alpha_delta_result) :: ad
   type(qq_result) :: qq
   type(ghyp_model_type) :: marginal, standardized
   real(dp) :: sc(2,2), a(1,2), es, omega

   model=ghyp_ad(0.7_dp,1.8_dp,1.2_dp,[0.3_dp],[0.2_dp],reshape([1.0_dp],[1,1]))
   call assert_close(ghyp_moment(model,1.0_dp,central=.false.),0.5209139881916593_dp,3.0e-11_dp,'raw mean')
   call assert_close(ghyp_skewness(model),0.3710470595901663_dp,2.0e-9_dp,'skewness')
   call assert_close(ghyp_kurtosis(model),4.371328532133212_dp,2.0e-8_dp,'kurtosis')
   es=esghyp(0.95_dp,model,loss=.true.)
   call assert_close(es,3.008893693572547_dp,2.0e-6_dp,'expected shortfall')
   omega=ghyp_omega(0.0_dp,model)
   call assert_close(omega,3.776544671407474_dp,3.0e-6_dp,'omega')

   sc=reshape([1.0_dp,0.3_dp,0.3_dp,0.8_dp],[2,2])
   mv=ghyp_mv(0.7_dp,1.4_dp,2.3_dp,[0.1_dp,-0.2_dp],sc,[0.2_dp,-0.1_dp])
   a=reshape([0.4_dp,0.6_dp],[1,2])
   transformed=transform_ghyp(mv,a,[0.1_dp])
   call assert_true(transformed%ok .and. transformed%dimension()==1,'linear transform')
   call assert_close(transformed%mu(1),0.02_dp,2.0e-13_dp,'transformed location')
   marginal=subset_ghyp(mv,[2])
   call assert_true(marginal%ok .and. marginal%dimension()==1,'marginal extraction')
   standardized=standardize_ghyp(model)
   call assert_true(standardized%ok,'standardization')
   ad=ghyp_alpha_delta(model)
   call assert_true(ad%ok,'alpha-delta conversion')
   call assert_close(ad%alpha,1.8_dp,2.0e-10_dp,'alpha-delta alpha')
   qq=qqghyp_data([0.2_dp,-0.1_dp,0.7_dp],model)
   call assert_true(qq%ok .and. abs(qq%sample(1)+0.1_dp)<1.0e-14_dp,'QQ data')
   attrib=esghyp_attribution(0.9_dp,mv,[0.4_dp,0.6_dp],loss=.true.)
   call assert_true(attrib%ok .and. size(attrib%contribution)==2,'ES attribution')
   print '(a)', 'test_risk_transform: PASS'
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
end program test_risk_transform
