program test_em
   use r_compat, only: dp
   use matrixdist_em
   use matrixdist_bivariate_fit
   use matrixdist_multi_fit
   use matrixdist_mphstar_fit
   implicit none
   real(dp)::a(1),s(1,1),obs(4),w(4),rate,alpha2(1),s11(1,1),s12(1,1),s22(1,1)
   real(dp)::bobs(4,2),s3(1,1,2),y(4,2)
   integer::dobs(4),bdobs(4,2),mdobs(4,2)
   logical::delta(4,2)
   real(dp)::rew(1,2),ys(4,2),tt(4)
   a=1.0_dp
   s=-1.0_dp
   obs=[0.5_dp,1.0_dp,2.0_dp,1.5_dp]
   w=1.0_dp
   call emstep_ph(a,s,obs,w)
   rate=4.0_dp/sum(obs)
   call chk(a(1),1.0_dp,1e-10_dp,'ph em alpha')
   call chk(s(1,1),-rate,2e-9_dp,'ph em rate')

   a=1.0_dp
   s=0.5_dp
   dobs=[1,2,4,3]
   call emstep_dph(a,s,dobs)
   call chk(s(1,1),real(sum(dobs)-size(dobs),dp)/real(sum(dobs),dp),2e-10_dp,'dph em q')

   a=1.0_dp
   s11=-1.0_dp
   s12=1.0_dp
   s22=-1.0_dp
   bobs(:,1)=obs
   bobs(:,2)=[0.3_dp,0.7_dp,1.1_dp,1.9_dp]
   call emstep_bivph(a,s11,s12,s22,bobs,w)
   call chk(s11(1,1),-4.0_dp/sum(bobs(:,1)),3e-8_dp,'bivph em s11')
   call chk(s12(1,1),4.0_dp/sum(bobs(:,1)),3e-8_dp,'bivph em s12')
   call chk(s22(1,1),-4.0_dp/sum(bobs(:,2)),3e-8_dp,'bivph em s22')

   a=1.0_dp
   s11=0.5_dp
   s12=0.5_dp
   s22=0.5_dp
   bdobs(:,1)=[1,2,4,3]
   bdobs(:,2)=[2,1,3,5]
   call emstep_bivdph(a,s11,s12,s22,bdobs,w)
   call chk(s11(1,1),real(sum(bdobs(:,1))-4,dp)/real(sum(bdobs(:,1)),dp),3e-9_dp,'bivdph em q1')
   call chk(s12(1,1),4.0_dp/real(sum(bdobs(:,1)),dp),3e-9_dp,'bivdph em cross')
   call chk(s22(1,1),real(sum(bdobs(:,2))-4,dp)/real(sum(bdobs(:,2)),dp),3e-9_dp,'bivdph em q2')

   alpha2=1.0_dp
   s3(1,1,:)=[-1.0_dp,-1.0_dp]
   y=bobs
   delta=.true.
   call emstep_mph_rc(alpha2,s3,y,delta,w)
   call chk(s3(1,1,1),-4.0_dp/sum(y(:,1)),4e-8_dp,'mph em rate1')
   call chk(s3(1,1,2),-4.0_dp/sum(y(:,2)),4e-8_dp,'mph em rate2')

   alpha2=1.0_dp
   s3(1,1,:)=[0.5_dp,0.5_dp]
   mdobs=bdobs
   call emstep_mdph(alpha2,s3,mdobs,w)
   call chk(s3(1,1,1),real(sum(mdobs(:,1))-4,dp)/real(sum(mdobs(:,1)),dp),4e-8_dp,'mdph em q1')
   call chk(s3(1,1,2),real(sum(mdobs(:,2))-4,dp)/real(sum(mdobs(:,2)),dp),4e-8_dp,'mdph em q2')

   alpha2=1.0_dp
   s1_mphstar: block
      real(dp)::ss(1,1)
      ss=-1.0_dp
      rew(1,:)=[0.3_dp,0.7_dp]
      tt=[0.5_dp,1.0_dp,2.0_dp,1.5_dp]
      ys(:,1)=0.3_dp*tt
      ys(:,2)=0.7_dp*tt
      call emstep_mphstar(alpha2,ss,rew,ys,reward_tol=1e-12_dp)
      call chk(ss(1,1),-4.0_dp/sum(tt),4e-8_dp,'mphstar em rate')
      call chk(rew(1,1),0.3_dp,4e-8_dp,'mphstar reward1')
      call chk(rew(1,2),0.7_dp,4e-8_dp,'mphstar reward2')
   end block s1_mphstar
   print *, 'test_em: PASS'
contains
   subroutine chk(got,want,eps,label)
      real(dp),intent(in)::got,want,eps
      character(len=*),intent(in)::label
      if(abs(got-want)>eps*max(1.0_dp,abs(want))) then
      print *,trim(label),got,want
      error stop 1
      end if
   end subroutine
end program
