program test_l1pack
  use l1pack
  implicit none
  integer :: fails, info, i
  real(dp), parameter :: tol=5.0e-6_dp
  real(dp) :: x1(5,1), y1(5), x2(5,2), p, q, d, m, expected
  real(dp) :: center(2), scatter(2,2), samp(4000,2), sm(2), cv(2,2)
  real(dp) :: data(8,2), z(8), ci(2,2), simlad(5,2)
  type(l1fit_result) :: l1
  type(lad_result) :: br, em
  type(laplace_fit_result) :: lf
  type(spatial_median_result) :: sf
  type(l1ccc_result) :: cc
  type(laplace_envelope_result) :: env

  fails=0
  call check_close(dlaplace(0.0_dp),1.0_dp/sqrt(2.0_dp),1e-12_dp,'dlaplace(0)',fails)
  call check_close(plaplace(0.0_dp),0.5_dp,1e-12_dp,'plaplace(0)',fails)
  do i=1,9
    p=real(i,dp)/10.0_dp
    q=qlaplace(p,location=1.2_dp,scale=2.3_dp)
    call check_close(plaplace(q,location=1.2_dp,scale=2.3_dp),p,2e-13_dp,'Laplace inversion',fails)
  end do

  x1(:,1)=[0.0_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp]
  y1=[1.0_dp,3.0_dp,5.0_dp,7.0_dp,100.0_dp]
  call l1fit(x1,y1,l1,intercept=.true.)
  call check_close(l1%coefficients(1),1.0_dp,tol,'l1fit intercept',fails)
  call check_close(l1%coefficients(2),2.0_dp,tol,'l1fit slope',fails)
  call check_close(l1%minimum,91.0_dp,tol,'l1fit SAD',fails)

  x2(:,1)=1.0_dp;x2(:,2)=x1(:,1)
  call lad_fit_br(x2,y1,br)
  call check_close(br%coefficients(1),1.0_dp,tol,'BR intercept',fails)
  call check_close(br%coefficients(2),2.0_dp,tol,'BR slope',fails)
  call check_close(br%sad,91.0_dp,tol,'BR SAD',fails)
  if(br%scale<=0.0_dp.or.br%loglik>=0.0_dp)call fail('BR scale/logLik',fails)
  call confint_lad(br,0.95_dp,ci)
  if(any(ci(:,1)>br%coefficients).or.any(ci(:,2)<br%coefficients))call fail('BR confidence intervals',fails)
  call simulate_lad(br,2,simlad)
  if(any(abs(simlad)>huge(1.0_dp)))call fail('simulate_lad finite',fails)

  call lad_fit_em(x2,y1,em,tol=1.0e-8_dp,maxiter=500)
  if(abs(em%coefficients(1)-1.0_dp)>5e-3_dp.or.abs(em%coefficients(2)-2.0_dp)>5e-3_dp) &
    call fail('EM coefficients',fails)
  if(em%sad>91.01_dp)call fail('EM objective',fails)

  center=[0.5_dp,-1.0_dp]
  scatter=reshape([1.0_dp,0.25_dp,0.25_dp,2.0_dp],[2,2])
  call rmlaplace(size(samp,1),center,scatter,samp,info)
  if(info/=0)call fail('rmLaplace info',fails)
  sm=sum(samp,dim=1)/real(size(samp,1),dp)
  if(maxval(abs(sm-center))>0.35_dp)call fail('rmLaplace mean',fails)
  cv=0.0_dp
  do i=1,size(samp,1)
    cv(1,:)=cv(1,:)+(samp(i,1)-sm(1))*(samp(i,:)-sm)
    cv(2,:)=cv(2,:)+(samp(i,2)-sm(2))*(samp(i,:)-sm)
  end do
  cv=cv/real(size(samp,1),dp)
  ! For this parameterization Cov(X)=4(p+1) Scatter = 12 Scatter for p=2.
  if(maxval(abs(cv/12.0_dp-scatter))>0.35_dp)call fail('rmLaplace covariance',fails)
  d=dmlaplace(center,center,scatter)
  if(d<=0.0_dp)call fail('dmLaplace positive',fails)

  data(:,1)=[-2.0_dp,-1.0_dp,0.0_dp,0.5_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp]
  data(:,2)=[-1.5_dp,-0.5_dp,0.2_dp,0.8_dp,1.2_dp,2.5_dp,2.7_dp,4.3_dp]
  call laplace_fit(data,lf,tol=1.0e-7_dp,maxiter=100)
  if(.not.lf%converged)call fail('LaplaceFit convergence',fails)
  if(any(lf%weights<=0.0_dp))call fail('LaplaceFit weights',fails)
  expected=0.0_dp
  do i=1,size(data,1)
    expected=expected+log_dmlaplace(data(i,:),lf%center,lf%scatter)
  end do
  call check_close(lf%loglik,expected,1e-9_dp,'LaplaceFit logLik identity',fails)
  z=wh_laplace(lf%distances,2)
  if(any(abs(z)>huge(1.0_dp)))call fail('WH finite',fails)

  call spatial_median_fit(data,sf,tol=1.0e-6_dp,maxiter=100)
  if(.not.sf%converged)call fail('spatial median convergence',fails)
  if(any(sf%weights<=0.0_dp))call fail('spatial median weights',fails)

  call l1ccc(data,cc,equal_means=.false.,boots=.false.)
  if(abs(cc%rho1)>1.5_dp.or.abs(cc%gaussian_rho1)>1.5_dp)call fail('CCC range',fails)
  call ustat_rho1(data,cc%ustat)
  if(cc%ustat%var_rho1<0.0_dp)call fail('U-stat variance',fails)

  ! Perfect agreement must produce rho1=1 for all three definitions.
  data(:,2)=data(:,1)
  call l1ccc(data,cc,equal_means=.false.,boots=.false.)
  call check_close(cc%rho1,1.0_dp,2e-5_dp,'Laplace rho1 perfect agreement',fails)
  call check_close(cc%gaussian_rho1,1.0_dp,2e-12_dp,'Gaussian rho1 perfect agreement',fails)
  call check_close(cc%ustat%rho1,1.0_dp,2e-12_dp,'Ustat rho1 perfect agreement',fails)

  ! Small smoke test for the non-plotting simulated envelope.
  data(:,2)=[-1.5_dp,-0.5_dp,0.2_dp,0.8_dp,1.2_dp,2.5_dp,2.7_dp,4.3_dp]
  call laplace_fit(data,lf,tol=1e-5_dp,maxiter=50)
  call envelope_laplace(data,lf,reps=4,conf=0.8_dp,out=env)
  if(size(env%lower)/=size(data,1).or.any(env%lower>env%upper))call fail('Laplace envelope',fails)

  ! Basic checks for nuisance estimate and LAD quantile residual transformation.
  z(1:5)=lad_quantile_residuals(br,y1)
  if(any(abs(z(1:5))>huge(1.0_dp)))call fail('quantile residual finite',fails)
  m=nuisance_vcov([ -3.0_dp,-1.0_dp,0.5_dp,2.0_dp,4.0_dp ],0.05_dp,.true.)
  if(m<=0.0_dp)call fail('nuisance_vcov positive',fails)

  if(fails==0)then
    print '(a)', 'test_l1pack: PASS'
  else
    print '(a,i0)', 'test_l1pack: FAIL ',fails
    error stop 1
  end if
contains
  subroutine check_close(a,b,t,name,fails)
    real(dp),intent(in)::a,b,t
    character(len=*),intent(in)::name
    integer,intent(inout)::fails
    if(abs(a-b)>t)then
      print '(a,2es20.10)',trim(name)//' mismatch: ',a,b
      fails=fails+1
    end if
  end subroutine
  subroutine fail(name,fails)
    character(len=*),intent(in)::name
    integer,intent(inout)::fails
    print '(a)',trim(name)//': FAIL'
    fails=fails+1
  end subroutine
end program test_l1pack
