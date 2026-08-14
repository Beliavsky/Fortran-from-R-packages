program test_density_sim
  use mclust
  implicit none
  integer,parameter::n=100
  real(dp)::x(n,1),p(3),q(3),cdf(3),lev(2)
  real(dp),allocatable::ld(:),xs(:,:)
  integer,allocatable::comp(:)
  type(mclust_fit)::fit
  integer::i,st
  do i=1,n/2; x(i,1)=-2.0_dp+0.25_dp*sin(0.73_dp*i); end do
  do i=n/2+1,n; x(i,1)=2.5_dp+0.35_dp*cos(0.47_dp*i); end do
  call fit_model(x,2,'V',fit)
  if(fit%status<0) error stop '1d fit'
  p=[0.1_dp,0.5_dp,0.9_dp]
  call quantile_mclust_1d(fit,p,q,st); if(st/=0)error stop 'quantile'
  call cdf_mclust_1d(fit,q,cdf,st); if(st/=0)error stop 'cdf'
  if(maxval(abs(cdf-p))>1e-7_dp) then; print *,cdf,p,q; error stop 'cdf-quantile mismatch'; end if
  call mixture_log_density(fit,x,ld,st); if(st/=0 .or. size(ld)/=n)error stop 'density'
  call hdr_levels(exp(ld),[0.5_dp,0.9_dp],lev)
  if(lev(1)<lev(2)) error stop 'hdr order'
  call simulate_fit(fit,500,xs,comp,st); if(st/=0)error stop 'simulation'
  if(any(comp<1) .or. any(comp>2))error stop 'component labels'
  if(abs(sum(real([(count(comp==i),i=1,2)],dp))/500.0_dp-1.0_dp)>1e-12_dp)error stop 'component counts'
  print *, 'test_density_sim PASS ', q
end program test_density_sim
