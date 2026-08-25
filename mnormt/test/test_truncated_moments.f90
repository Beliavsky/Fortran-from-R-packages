program test_truncated_moments
  use mnormt
  implicit none
  real(dp) :: mu1(1),s1(1,1),lo1(1),up1(1),x1(1),p,d,vref
  real(dp) :: mu2(2),s2(2,2),lo2(2),up2(2),samp(3000,2),sm(2)
  integer :: kap2(2),info,fails
  real(dp), allocatable :: raw(:)
  type(trunc_moment_result) :: mr
  fails=0
  mu1=0.0_dp;s1(1,1)=1.0_dp;lo1=-1.0_dp;up1=1.0_dp;x1=0.0_dp
  d=dmtruncnorm(x1,mu1,s1,lo1,up1)
  p=normal_pdf(0.0_dp)/(normal_cdf(1.0_dp)-normal_cdf(-1.0_dp))
  if(abs(d-p)>1e-12_dp) then;print*,'dmtrunc1 fail',d,p;fails=fails+1;end if
  p=pmtruncnorm([0.0_dp],mu1,s1,lo1,up1)
  if(abs(p-0.5_dp)>1e-12_dp) then;print*,'pmtrunc1 fail',p;fails=fails+1;end if
  mu2=0.0_dp;s2=0.0_dp;s2(1,1)=1.0_dp;s2(2,2)=1.0_dp;lo2=[-1.0_dp,-2.0_dp];up2=[1.0_dp,2.0_dp];kap2=[2,2]
  call recintab(kap2,lo2,up2,mu2,s2,raw,info)
  if(info/=0) then;print*,'recintab status',info;fails=fails+1;end if
  mr=mom_mtruncnorm(kap2,mu2,s2,lo2,up2)
  vref=1.0_dp-2.0_dp*normal_pdf(1.0_dp)/(normal_cdf(1.0_dp)-normal_cdf(-1.0_dp))
  if(maxval(abs(mr%mean))>1e-12_dp) then;print*,'moment mean fail',mr%mean;fails=fails+1;end if
  if(abs(mr%covariance(1,1)-vref)>2e-10_dp) then;print*,'moment var fail',mr%covariance(1,1),vref;fails=fails+1;end if
  if(abs(mr%covariance(1,2))>1e-12_dp) then;print*,'moment cov fail',mr%covariance(1,2);fails=fails+1;end if
  call rmtruncnorm(3000,mu2,s2,lo2,up2,samp,burnin=20,info=info)
  sm=sum(samp,dim=1)/3000.0_dp
  if(maxval(abs(sm))>0.08_dp .or. any(samp<=spread(lo2,1,3000)) .or. any(samp>=spread(up2,1,3000))) then
    print*,'rmtrunc fail',sm;fails=fails+1
  end if
  ! Univariate t specialization and truncation.
  p=pmt([0.5_dp],[0.0_dp],reshape([1.0_dp],[1,1]),5.0_dp)
  if(abs(p-0.6808505642_dp)>2e-9_dp) then;print*,'pmt1 fail',p;fails=fails+1;end if
  p=pmtrunct([0.0_dp],[0.0_dp],reshape([1.0_dp],[1,1]),5.0_dp,[-1.0_dp],[1.0_dp])
  if(abs(p-0.5_dp)>1e-10_dp) then;print*,'pmtrunct fail',p;fails=fails+1;end if
  if(fails==0) then; print*,'test_truncated_moments: PASS'; else; print*,'test_truncated_moments: FAIL',fails; error stop 1; end if
end program
