program test_fit
  use dirichletreg, only : dp, design_block, dirichletreg_model, seed_rng, rdirichlet, fit_common, fit_alternative
  implicit none
  integer, parameter :: n=1200,d=3
  real(dp) :: alpha0(d), y(n,d), x(n,1), z(n,1), true_alt(3), mu0(d), phi0
  type(design_block) :: xb(d)
  type(dirichletreg_model) :: mc,ma
  integer :: j,stat,failures

  failures=0
  alpha0=[2.0_dp,3.0_dp,5.0_dp]
  call seed_rng(777)
  call rdirichlet(n,alpha0,y,stat)
  do j=1,d
    allocate(xb(j)%x(n,1)); xb(j)%x(:,1)=1.0_dp
  end do
  call fit_common(y,xb,mc,iterlim=1000)
  if(mc%convergence>2) failures=failures+1
  if(maxval(abs(exp(mc%coefficients)-alpha0))>0.45_dp) failures=failures+1

  mu0=[0.2_dp,0.3_dp,0.5_dp]; phi0=12.0_dp; alpha0=mu0*phi0
  call seed_rng(778); call rdirichlet(n,alpha0,y,stat)
  x(:,1)=1.0_dp; z(:,1)=1.0_dp
  true_alt=[log(mu0(1)/mu0(3)),log(mu0(2)/mu0(3)),log(phi0)]
  call fit_alternative(y,x,z,3,ma,iterlim=1000)
  if(ma%convergence>2) failures=failures+1
  if(maxval(abs(ma%coefficients-true_alt))>0.20_dp) failures=failures+1

  if(failures==0) then
    print '(a)', 'test_fit: PASS'
  else
    print '(a,i0)', 'test_fit: FAIL ',failures
    print *, 'common coeff=',mc%coefficients,' conv=',mc%convergence
    print *, 'alt coeff=',ma%coefficients,' true=',true_alt,' conv=',ma%convergence
    error stop 1
  end if
end program test_fit
