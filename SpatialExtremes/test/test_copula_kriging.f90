program test_copula_kriging
   use SpatialExtremes
   implicit none
   real(dp)::u(3,2),cor(2,2),ll,coord(2,1),newc(1,1),w(2,1),obs(2),pred(1)
   integer::info
   u=reshape([0.2_dp,0.4_dp,0.7_dp,0.8_dp,0.6_dp,0.3_dp],[3,2])
   cor=0.0_dp
   cor(1,1)=1
   cor(2,2)=1
   ll=gaussian_copula_loglik(u,cor)
   call check(abs(ll)<1e-11_dp,'Gaussian independence copula')
   ll=student_copula_loglik(u,5.0_dp,cor)
   call check(abs(ll-0.15382289912120184_dp)<2e-10_dp,'Student copula reference')
   coord(:,1)=[0.0_dp,2.0_dp]
   newc(1,1)=1.0_dp
   obs=[1.0_dp,3.0_dp]
   call simple_kriging_weights(coord,newc,COV_POWEREXP,1.0_dp,1.0_dp,1.0_dp,0.0_dp,w,info)
   call check(info==0,'kriging solve')
   pred=simple_kriging_predict(obs,w)
   call check(abs(w(1,1)-w(2,1))<1e-13_dp,'symmetric kriging weights')
   call check(pred(1)>0.0_dp,'kriging prediction')
   print *,'test_copula_kriging: PASS'
contains
   subroutine check(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then
      print *,'FAIL: ',msg
      error stop 1
      end if
   end subroutine
end program
