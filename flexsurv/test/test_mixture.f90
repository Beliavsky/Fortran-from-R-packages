program test_mixture
  use flexsurv
  implicit none
  type(flexsurvmix_result)::m,fit
  type(flexsurv_spec)::sp(2)
  type(flexsurv_data)::dat
  real(dp)::s,f,exact,t(20)
  integer::st(20),i,fails
  fails=0
  call initialize_spec(sp(1),dist_exponential,1,[2.0_dp])
  call initialize_spec(sp(2),dist_exponential,1,[0.3_dp])
  allocate(m%mixing(2),m%specs(2),m%components(2));m%mixing=[0.4_dp,0.6_dp];m%specs=sp
  allocate(m%components(1)%theta(1),m%components(2)%theta(1))
  m%components(1)%theta=log(2.0_dp);m%components(2)%theta=log(0.3_dp)
  s=mix_survival(m,1,1.0_dp);exact=0.4_dp*exp(-2.0_dp)+0.6_dp*exp(-0.3_dp)
  if(abs(s-exact)>1e-12_dp)fails=fails+1
  f=mix_density(m,1,1.0_dp);exact=0.4_dp*2.0_dp*exp(-2.0_dp)+0.6_dp*0.3_dp*exp(-0.3_dp)
  if(abs(f-exact)>1e-12_dp)fails=fails+1
  ! Synthetic separated sample, deterministic exponential quantiles.
  do i=1,10;t(i)=-log(1.0_dp-(real(i,dp)-0.5_dp)/10.0_dp)/2.0_dp;end do
  do i=11,20;t(i)=-log(1.0_dp-(real(i-10,dp)-0.5_dp)/10.0_dp)/0.3_dp;end do
  st=1;call prepare_survival_data(dat,t,st)
  call initialize_spec(sp(1),dist_exponential,20,[1.8_dp])
  call initialize_spec(sp(2),dist_exponential,20,[0.35_dp])
  fit=fit_flexsurvmix(dat,sp,[0.5_dp,0.5_dp],maxit=60,tol=1e-7_dp)
  if(.not.ieee_finite(fit%loglik))fails=fails+1
  if(abs(sum(fit%mixing)-1.0_dp)>1e-12_dp)fails=fails+1
  if(any(fit%mixing<=0.0_dp))fails=fails+1
  if(fails>0)then;print *,'mix ',fit%mixing,fit%loglik;error stop 1;end if
  print *,'test_mixture: PASS'
contains
  logical function ieee_finite(x)
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    real(dp),intent(in)::x
    ieee_finite=ieee_is_finite(x)
  end function ieee_finite
end program test_mixture
