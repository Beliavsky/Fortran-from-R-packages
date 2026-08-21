program test_spectral_methods
  use tsa, only : dp, spectral_estimate, spec_ar, spec_pgram, &
    tskernel_weights, daniell_kernel, modified_daniell_kernel, &
    fejer_kernel, dirichlet_kernel, modified_daniell_weights
  implicit none
  integer :: failures, i, st, st2
  real(dp) :: e, u
  real(dp), allocatable :: x(:), xm(:,:), compact(:), full(:), full2(:), k(:)
  type(spectral_estimate) :: syw, sb, sb2, so, sm, sa, p1, p2, pm

  failures=0
  allocate(x(500))
  x=0.0_dp
  do i=3,size(x)
    u=sin(12.9898_dp*real(i,dp))*43758.5453123_dp
    e=(u-floor(u))-0.5_dp
    x(i)=0.65_dp*x(i-1)-0.20_dp*x(i-2)+e
  end do

  syw=spec_ar(x,n_freq=33,order=2,method='yule-walker')
  sb =spec_ar(x,n_freq=33,order=2,method='burg',var_method=1)
  sb2=spec_ar(x,n_freq=33,order=2,method='burg',var_method=2)
  so =spec_ar(x,n_freq=33,order=2,method='ols')
  sm =spec_ar(x,n_freq=33,order=2,method='mle')
  call check(syw%status==0.and.sb%status==0.and.sb2%status==0.and.so%status==0.and.sm%status==0, &
    'all spec.ar methods fit fixed AR(2)',failures)
  call check(all(syw%spectrum>0.0_dp).and.all(sb%spectrum>0.0_dp).and.all(so%spectrum>0.0_dp).and. &
    all(sm%spectrum>0.0_dp),'all AR spectra positive',failures)
  call check(maxval(abs(so%spectrum-sm%spectrum))/maxval(sm%spectrum)<0.20_dp, &
    'OLS and MLE spectra agree on long AR sample',failures)
  call check(maxval(abs(sb%spectrum-sb2%spectrum))/maxval(sb%spectrum)<0.20_dp, &
    'Burg variance methods retain common AR shape',failures)

  sa=spec_ar(x,n_freq=17,order_max=6,method='burg')
  call check(sa%status==0.and.sa%order>=1.and.sa%order<=6,'Burg AIC order selection',failures)
  sa=spec_ar(x,n_freq=17,order_max=6,method='ols')
  call check(sa%status==0.and.sa%order>=1.and.sa%order<=6,'OLS AIC order selection',failures)
  sa=spec_ar(x,n_freq=17,order_max=4,method='mle')
  call check(sa%status==0.and.sa%order>=1.and.sa%order<=4,'MLE AIC order selection',failures)

  compact=modified_daniell_kernel(2)
  call tskernel_weights(compact,full,st)
  call check(st==0.and.size(full)==5,'compact tskernel expansion',failures)
  call check(maxval(abs(full-[0.125_dp,0.25_dp,0.25_dp,0.25_dp,0.125_dp]))<1.0e-14_dp, &
    'modified Daniell coefficients match R convention',failures)
  p1=spec_pgram(x,kernel_coef=compact,taper=0.1_dp,fast=.false.)
  p2=spec_pgram(x,kernel_weights=full,taper=0.1_dp,fast=.false.)
  call check(p1%status==0.and.p2%status==0.and.maxval(abs(p1%spectrum-p2%spectrum))<1.0e-12_dp, &
    'compact and expanded kernels are equivalent',failures)

  k=daniell_kernel(3); call tskernel_weights(k,full,st)
  call check(st==0.and.abs(sum(full)-1.0_dp)<1.0e-14_dp,'Daniell kernel normalization',failures)
  k=fejer_kernel(4,3.5_dp); call tskernel_weights(k,full,st)
  call check(st==0.and.abs(sum(full)-1.0_dp)<1.0e-12_dp,'Fejer kernel normalization',failures)
  k=dirichlet_kernel(4,2.25_dp); call tskernel_weights(k,full,st)
  call check(st==0.and.abs(sum(full)-1.0_dp)<1.0e-12_dp,'Dirichlet kernel normalization',failures)

  k=modified_daniell_kernel([1,2])
  call tskernel_weights(k,full,st)
  call modified_daniell_weights([3,5],full2,st2)
  call check(st==0.and.st2==0.and.size(full)==size(full2), &
    'composite modified Daniell size',failures)
  if(st==0.and.st2==0.and.size(full)==size(full2))then
    call check(maxval(abs(full-full2))<1.0e-13_dp,'composite modified Daniell',failures)
  end if

  k=daniell_kernel([1,2])
  call tskernel_weights(k,full,st)
  call check(st==0.and.abs(sum(full)-1.0_dp)<1.0e-12_dp,'composite Daniell',failures)

  allocate(xm(size(x),2)); xm(:,1)=x; xm(:,2)=0.6_dp*x
  do i=2,size(x); xm(i,2)=xm(i,2)+0.2_dp*x(i-1); end do
  compact=modified_daniell_kernel(2)
  pm=spec_pgram(xm,kernel_coef=compact,taper_series=[0.0_dp,0.2_dp],fast=.false.)
  p1=spec_pgram(xm(:,1),kernel_coef=compact,taper=0.0_dp,fast=.false.)
  p2=spec_pgram(xm(:,2),kernel_coef=compact,taper=0.2_dp,fast=.false.)
  call check(pm%status==0.and.allocated(pm%taper_series).and.size(pm%taper_series)==2, &
    'per-series taper accepted',failures)
  call check(maxval(abs(pm%spectrum_matrix(:,1)-p1%spectrum))<1.0e-11_dp.and. &
    maxval(abs(pm%spectrum_matrix(:,2)-p2%spectrum))<1.0e-11_dp, &
    'multivariate taper autospectra equal separate univariate fits',failures)
  call check(allocated(pm%degrees_freedom_series).and. &
    abs(pm%degrees_freedom_series(1)-pm%degrees_freedom_series(2))>1.0e-6_dp, &
    'per-series taper degrees of freedom retained',failures)
  call check(all(pm%coherence>=-1.0e-12_dp).and.all(pm%coherence<=1.0_dp+1.0e-8_dp), &
    'coherence remains in unit interval',failures)

  if(failures==0)then
    print '(a)','test_spectral_methods: PASS'
  else
    print '(a,i0)','test_spectral_methods: FAIL ',failures
    error stop 1
  end if
contains
  subroutine check(ok,label,failures)
    logical,intent(in)::ok
    character(len=*),intent(in)::label
    integer,intent(inout)::failures
    if(.not.ok)then
      failures=failures+1
      print '(a,a)','FAIL: ',trim(label)
    end if
  end subroutine check
end program test_spectral_methods
