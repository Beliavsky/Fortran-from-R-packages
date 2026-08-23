program test_normalp
  use normalp
  implicit none
  integer, parameter :: n=6000
  real(dp), parameter :: pi=acos(-1.0_dp)
  real(dp) :: p0, q0, stats(3), err
  real(dp), allocatable :: z(:), xx(:,:), yy(:)
  integer, allocatable :: seed(:)
  type(normalp_params) :: fit
  type(lmp_result) :: reg
  integer :: i, ns

  call random_seed(size=ns); allocate(seed(ns)); seed=13579; call random_seed(put=seed)

  call check(abs(dnormp(0.0_dp,p=2.0_dp)-1.0_dp/sqrt(2.0_dp*pi))<1.0e-12_dp,'normal density')
  call check(abs(pnormp(0.0_dp,p=2.0_dp)-0.5_dp)<1.0e-14_dp,'median cdf')
  call check(abs(dnormp(0.0_dp,p=1.0_dp)-0.5_dp)<1.0e-12_dp,'laplace density')
  call check(abs(pnormp(log(2.0_dp),p=1.0_dp)-0.75_dp)<2.0e-11_dp,'laplace cdf')

  do i=1,9
    p0=real(i,dp)/10.0_dp
    q0=qnormp(p0,mu=1.2_dp,sigmap=0.7_dp,p=3.5_dp)
    err=abs(pnormp(q0,mu=1.2_dp,sigmap=0.7_dp,p=3.5_dp)-p0)
    call check(err<2.0e-11_dp,'cdf quantile inversion')
  end do

  allocate(z(n)); call rnormp(z,mu=2.0_dp,sigmap=1.5_dp,p=2.0_dp)
  call check(abs(sum(z)/real(n,dp)-2.0_dp)<0.08_dp,'rng mean')
  call check(abs(sqrt(sum((z-sum(z)/real(n,dp))**2)/real(n,dp))-1.5_dp)<0.08_dp,'rng sd p2')

  call rnormp(z,mu=-0.3_dp,sigmap=1.1_dp,p=3.0_dp)
  call paramp_fit(z,fit,p_fixed=3.0_dp)
  call check(abs(fit%mp+0.3_dp)<0.08_dp,'paramp location fixed p')
  call check(abs(fit%sp-1.1_dp)<0.08_dp,'paramp scale fixed p')
  p0=estimatep(z,-0.3_dp,2.0_dp)
  call check(p0>2.0_dp .and. p0<4.5_dp,'estimatep shape')

  call kurtosis_p(stats,p=2.0_dp)
  call check(abs(stats(2)-3.0_dp)<1.0e-12_dp,'normal kurtosis')
  call check(abs(stats(3)-3.0_dp)<1.0e-12_dp,'Bp identity')

  allocate(xx(300,2),yy(300));
  do i=1,300
    xx(i,1)=1.0_dp; xx(i,2)=-2.0_dp+4.0_dp*real(i-1,dp)/299.0_dp
  end do
  call rnormp(yy,mu=0.0_dp,sigmap=0.15_dp,p=2.0_dp)
  yy=1.25_dp-0.8_dp*xx(:,2)+yy
  call lmp_fit(xx,yy,reg,p_fixed=2.0_dp)
  call check(abs(reg%coef(1)-1.25_dp)<0.05_dp,'lmp intercept')
  call check(abs(reg%coef(2)+0.8_dp)<0.05_dp,'lmp slope')

  print '(a)', 'test_normalp: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok; character(*),intent(in)::msg
    if(.not.ok) then; print '(a,a)','FAIL: ',msg; error stop 1; end if
  end subroutine
end program
