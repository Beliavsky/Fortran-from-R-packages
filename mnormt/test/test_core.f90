program test_core
  use mnormt
  implicit none
  real(dp) :: mu2(2),s2(2,2),x2(2),p,dens,lo2(2),up2(2)
  real(dp) :: mu3(3),s3(3,3),x3(3),p3
  real(dp) :: sam(4000,2),m(2),v(2)
  type(probability_result) :: pr
  integer :: i,info,fails
  fails=0
  mu2=[0.0_dp,0.0_dp]; s2=reshape([1.0_dp,0.3_dp,0.3_dp,1.0_dp],[2,2]); x2=[0.5_dp,-0.25_dp]
  p=pmnorm(x2,mu2,s2)
  if(abs(p-0.3179461572_dp)>2e-8_dp) then; print *, 'pmnorm2 fail',p; fails=fails+1; end if
  lo2=[-huge(1.0_dp),-huge(1.0_dp)]; up2=x2
  pr=sadmvn_prob(lo2,up2,mu2,s2,20000,1e-10_dp,0.0_dp)
  if(abs(pr%value-p)>2e-6_dp) then; print *, 'sadmvn2 fail',pr%value,p; fails=fails+1; end if
  dens=dmnorm([0.0_dp,0.0_dp],mu2,s2)
  if(abs(dens-1.0_dp/(2.0_dp*acos(-1.0_dp)*sqrt(1.0_dp-0.3_dp**2)))>1e-12_dp) fails=fails+1
  p=biv_nt_prob(5.0_dp,lo2,up2,mu2,s2)
  if(abs(p-0.31574_dp)>5e-4_dp) then; print *, 'bivt fail',p; fails=fails+1; end if
  mu3=[0.0_dp,0.0_dp,0.0_dp]; s3=0.0_dp; do i=1,3;s3(i,i)=1.0_dp;end do
  x3=[0.0_dp,0.0_dp,0.0_dp]; p3=pmnorm(x3,mu3,s3)
  if(abs(p3-0.125_dp)>1e-10_dp) then; print *, 'pmnorm3 fail',p3; fails=fails+1; end if
  call rmnorm(4000,mu2,s2,sam,info)
  m=sum(sam,dim=1)/4000.0_dp
  v=[sum((sam(:,1)-m(1))**2),sum((sam(:,2)-m(2))**2)]/3999.0_dp
  if(maxval(abs(m))>0.08_dp .or. maxval(abs(v-1.0_dp))>0.08_dp) then; print *,'rng fail',m,v;fails=fails+1;end if
  if(fails==0) then; print *,'test_core: PASS'; else; print *,'test_core: FAIL',fails; error stop 1; end if
end program
