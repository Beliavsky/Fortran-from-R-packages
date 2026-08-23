program test_parity_v02
   use directional
   implicit none
   real(dp) :: ang(6), x(8,3), probs(2), kap(2), mus(2,3), ll
   real(dp) :: scmu(3), scrho, scll, pkmu(3), pkrho, pkll
   real(dp) :: gam(3), lam(3), lc(3), gmat(3,3), kd(2), a4(4,4), bx(20,4)
   real(dp) :: f(3,3), mf(3,3,8), u(3,3), s(3), v(3,3)
   integer :: pred(2), cl(8), i
   type(circular_mle_result) :: cm, ce
   type(vmf_mle_result) :: vm
   type(test_result) :: tr

   ang = [5.0_dp, 8.0_dp, 10.0_dp, 12.0_dp, 15.0_dp, 17.0_dp]
   cm = cardio_mle(ang)
   if (cm%rho < 0.0_dp .or. cm%rho > 0.5_dp) error stop 1
   ce = circexp_mle([0.1_dp,0.2_dp,0.3_dp,0.4_dp], rads=.true.)
   if (ce%lambda <= 0.0_dp) error stop 2

   x = reshape([ &
      0.99_dp,0.10_dp,0.02_dp, 0.98_dp,0.15_dp,0.05_dp, 0.97_dp,0.20_dp,0.04_dp, 0.96_dp,0.25_dp,0.03_dp, &
     -0.99_dp,-0.10_dp,0.02_dp,-0.98_dp,-0.15_dp,0.05_dp,-0.97_dp,-0.20_dp,0.04_dp,-0.96_dp,-0.25_dp,0.03_dp], [8,3], order=[2,1])
   do i=1,8; x(i,:)=x(i,:)/sqrt(sum(x(i,:)**2)); end do
   vm=vmf_mle(x(1:4,:)); if(vm%kappa<=0)error stop 3
   call mixvmf_mle(x,2,probs,kap,mus,ll,cl)
   if(abs(sum(probs)-1.0_dp)>1e-8_dp)error stop 4
   call dirknn(x([1,8],:),x,[1,1,1,1,2,2,2,2],3,pred)
   if(any(pred/=[1,2]))error stop 5
   call dirda_vmf(x([1,8],:),x,[1,1,1,1,2,2,2,2],pred)
   if(any(pred/=[1,2]))error stop 6

   tr=rayleigh_test(x(1:4,:),mc_reps=1); if(tr%statistic<=0 .or. tr%p_value<0 .or. tr%p_value>1)error stop 7
   tr=kuiper_test(ang); if(tr%statistic<=0 .or. tr%p_value<0 .or. tr%p_value>1)error stop 8

   call spcauchy_mle(x(1:4,:),scmu,scrho,scll,maxit=20)
   if(scrho<0 .or. scrho>=1)error stop 9
   call pkbd_mle(x(1:4,:),pkmu,pkrho,pkll,maxit=20)
   if(pkrho<0 .or. pkrho>=1)error stop 10

   gam=[0.0_dp,2.0_dp,0.0_dp];lam=[0.0_dp,-0.2_dp,0.2_dp];lc=fb_saddle(gam,lam)
   if(any(.not.(lc<huge(1.0_dp))))error stop 11
   gmat=0;gmat(1,1)=1;gmat(2,2)=1;gmat(3,3)=1
   kd=dkent(x(1:2,:),gmat,2.0_dp,0.2_dp)
   if(any(kd<=0))error stop 12

   a4=0;a4(2,2)=0.2_dp;a4(3,3)=0.5_dp;a4(4,4)=1.0_dp;bx=rbingham(20,a4)
   do i=1,20;if(abs(sum(bx(i,:)**2)-1)>1e-8_dp)error stop 13;end do
   f=0;f(1,1)=1.0_dp;f(2,2)=0.5_dp;f(3,3)=0.2_dp;mf=rmatrixfisher(8,f)
   call matrixfisher_mle(mf,u,s,v)
   if(any(s<0))error stop 14
   print *, 'test_parity_v02: PASS'
end program
