program test_derivatives
  use dirichletreg, only : dp, design_block, common_loglik_score, common_loglik_score_hessian, &
       alternative_loglik_score, alternative_loglik_score_hessian
  implicit none
  integer, parameter :: n=6,d=3,p=2,q=2
  real(dp) :: y(n,d), xv(n), w(n), theta_c(6), gc(6), hc(6,6), f, fp, fm, eps
  real(dp) :: gnc(6), hnc(6,6), tp(6), tm(6), gp(6), gm(6)
  real(dp) :: x(n,p), z(n,q), theta_a((d-1)*p+q), ga((d-1)*p+q), ha((d-1)*p+q,(d-1)*p+q)
  real(dp) :: gna((d-1)*p+q), hna((d-1)*p+q,(d-1)*p+q), ta((d-1)*p+q), tb((d-1)*p+q)
  type(design_block) :: xb(d)
  integer :: i,j,failures

  failures=0; eps=2.0e-6_dp
  xv=[-1.0_dp,-0.6_dp,-0.2_dp,0.2_dp,0.6_dp,1.0_dp]
  y=reshape([ &
    0.20_dp,0.25_dp,0.35_dp,0.30_dp,0.15_dp,0.40_dp, &
    0.30_dp,0.45_dp,0.25_dp,0.20_dp,0.50_dp,0.35_dp, &
    0.50_dp,0.30_dp,0.40_dp,0.50_dp,0.35_dp,0.25_dp],[n,d])
  w=[1.0_dp,0.8_dp,1.2_dp,1.0_dp,1.5_dp,0.7_dp]
  do j=1,d
    allocate(xb(j)%x(n,p)); xb(j)%x(:,1)=1.0_dp; xb(j)%x(:,2)=xv
  end do
  theta_c=[0.3_dp,0.1_dp,0.6_dp,-0.15_dp,0.9_dp,0.05_dp]
  call common_loglik_score_hessian(theta_c,y,xb,w,f,gc,hc)
  do j=1,size(theta_c)
    tp=theta_c; tm=theta_c; tp(j)=tp(j)+eps; tm(j)=tm(j)-eps
    call common_loglik_score(tp,y,xb,w,fp,gp)
    call common_loglik_score(tm,y,xb,w,fm,gm)
    gnc(j)=(fp-fm)/(2.0_dp*eps)
    hnc(:,j)=(gp-gm)/(2.0_dp*eps)
  end do
  if(maxval(abs(gc-gnc))>2.0e-6_dp) failures=failures+1
  if(maxval(abs(hc-hnc))>2.0e-5_dp) failures=failures+1

  x(:,1)=1.0_dp; x(:,2)=xv
  z(:,1)=1.0_dp; z(:,2)=xv*xv
  theta_a=[-0.4_dp,0.10_dp,-0.1_dp,-0.08_dp,1.4_dp,0.2_dp]
  call alternative_loglik_score_hessian(theta_a,y,x,z,3,w,f,ga,ha)
  do j=1,size(theta_a)
    ta=theta_a; tb=theta_a; ta(j)=ta(j)+eps; tb(j)=tb(j)-eps
    call alternative_loglik_score(ta,y,x,z,3,w,fp,gp)
    call alternative_loglik_score(tb,y,x,z,3,w,fm,gm)
    gna(j)=(fp-fm)/(2.0_dp*eps)
    hna(:,j)=(gp-gm)/(2.0_dp*eps)
  end do
  if(maxval(abs(ga-gna))>3.0e-6_dp) failures=failures+1
  if(maxval(abs(ha-hna))>5.0e-5_dp) then
    print '(a,es12.4)', 'alternative Hessian max error: ',maxval(abs(ha-hna))
    failures=failures+1
  end if

  if(failures==0) then
    print '(a)', 'test_derivatives: PASS'
  else
    print '(a,i0)', 'test_derivatives: FAIL ',failures
    error stop 1
  end if
end program test_derivatives
