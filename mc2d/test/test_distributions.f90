program test_distributions
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use mc2d, only : dp, pi, dbern, pbern, qbern, dbetagen, pbetagen, qbetagen, &
    dlnormb, plnormb, qlnormb, dpert, ppert, qpert, dpert_mean, dtriang, &
    ptriang, qtriang, dtriang_mean, dmqi, pmqi, qmqi, dempiricald, &
    pempiricald, qempiricald, dempiricalc, pempiricalc, qempiricalc, &
    ddirichlet, rdirichlet, dmultinomial, rmultinomial, dmultinormal, &
    rmultinormal
  implicit none
  integer :: fails, i
  real(dp) :: p, q, d, mn, expected
  real(dp), allocatable :: xmv(:,:)
  real(dp) :: xd(3), alpha(3), mqi(3), intrinsic(2), vals(3), probs(3)
  real(dp) :: mu2(2), sig2(2,2)
  integer :: counts(3)

  fails=0
  call check_close(dbern(0.0_dp,0.3_dp),0.7_dp,1e-14_dp,'dbern(0)',fails)
  call check_close(dbern(1.0_dp,0.3_dp),0.3_dp,1e-14_dp,'dbern(1)',fails)
  call check_close(pbern(0.0_dp,0.3_dp),0.7_dp,1e-14_dp,'pbern',fails)
  call check_close(qbern(0.7_dp,0.3_dp),0.0_dp,0.0_dp,'qbern boundary',fails)

  p=0.37_dp
  q=qbetagen(p,3.0_dp,5.0_dp,1.0_dp,6.0_dp)
  call check_close(pbetagen(q,3.0_dp,5.0_dp,1.0_dp,6.0_dp),p,2e-9_dp, &
    'betagen inverse',fails)
  call check_close(dbetagen(0.5_dp,2.0_dp,3.0_dp),1.5_dp,2e-12_dp, &
    'betagen beta identity',fails)

  expected=1.0_dp/sqrt(2.0_dp*pi)
  call check_close(dlnormb(1.0_dp),expected,2e-12_dp,'lognormalb defaults',fails)
  q=qlnormb(0.73_dp,5.0_dp,3.0_dp)
  call check_close(plnormb(q,5.0_dp,3.0_dp),0.73_dp,2e-9_dp, &
    'lognormalb inverse',fails)

  mn=(3.0_dp+4.0_dp*5.0_dp+10.0_dp)/6.0_dp
  call check_close(dpert(6.0_dp,3.0_dp,5.0_dp,10.0_dp), &
    dpert_mean(6.0_dp,3.0_dp,mn,10.0_dp),2e-12_dp,'PERT mean parametrization',fails)
  q=qpert(0.41_dp,3.0_dp,5.0_dp,10.0_dp)
  call check_close(ppert(q,3.0_dp,5.0_dp,10.0_dp),0.41_dp,2e-9_dp, &
    'PERT inverse',fails)
  if (.not.ieee_is_nan(dpert(2.0_dp,2.0_dp,2.0_dp,2.0_dp))) then
    print '(a)','FAIL PERT degenerate NaN'; fails=fails+1
  end if

  mn=(3.0_dp+6.0_dp+10.0_dp)/3.0_dp
  call check_close(dtriang(5.0_dp,3.0_dp,6.0_dp,10.0_dp), &
    dtriang_mean(5.0_dp,3.0_dp,mn,10.0_dp),2e-12_dp, &
    'triangular mean parametrization',fails)
  q=qtriang(0.61_dp,3.0_dp,6.0_dp,10.0_dp)
  call check_close(ptriang(q,3.0_dp,6.0_dp,10.0_dp),0.61_dp,2e-12_dp, &
    'triangular inverse',fails)

  mqi=[40.0_dp,50.0_dp,60.0_dp]; intrinsic=[0.0_dp,100.0_dp]
  q=qmqi(0.72_dp,mqi,intrinsic=intrinsic)
  call check_close(pmqi(q,mqi,intrinsic=intrinsic),0.72_dp,2e-12_dp, &
    'MQI inverse',fails)
  if (dmqi(50.0_dp,mqi,intrinsic=intrinsic)<0.0_dp) then
    print '(a)','FAIL MQI density'; fails=fails+1
  end if

  vals=[10.0_dp,20.0_dp,20.0_dp]; probs=[1.0_dp,2.0_dp,5.0_dp]
  call check_close(dempiricald(20.0_dp,vals,probs),7.0_dp/8.0_dp,1e-14_dp, &
    'empiricalD duplicate aggregation',fails)
  call check_close(pempiricald(10.0_dp,vals,probs),1.0_dp/8.0_dp,1e-14_dp, &
    'empiricalD CDF',fails)
  call check_close(qempiricald(0.5_dp,vals,probs),20.0_dp,0.0_dp, &
    'empiricalD quantile',fails)

  vals=[1.0_dp,3.0_dp,5.0_dp]; probs=[1.0_dp,2.0_dp,1.0_dp]
  do i=1,9
    p=real(i,dp)/10.0_dp
    q=qempiricalc(p,0.0_dp,6.0_dp,vals,probs)
    call check_close(pempiricalc(q,0.0_dp,6.0_dp,vals,probs),p,2e-10_dp, &
      'empiricalC inverse',fails)
  end do
  if (dempiricalc(-1.0_dp,0.0_dp,6.0_dp,vals,probs)/=0.0_dp) then
    print '(a)','FAIL empiricalC support'; fails=fails+1
  end if

  xd=[0.2_dp,0.3_dp,0.5_dp]; alpha=[1.0_dp,2.0_dp,3.0_dp]
  if (ddirichlet(xd,alpha)<0.0_dp) then
    print '(a)','FAIL Dirichlet density'; fails=fails+1
  end if
  call rdirichlet(alpha,xd)
  call check_close(sum(xd),1.0_dp,2e-14_dp,'Dirichlet sample sum',fails)

  counts=[1,2,3]
  d=dmultinomial(counts,[0.2_dp,0.3_dp,0.5_dp])
  expected=60.0_dp*0.2_dp*0.3_dp**2*0.5_dp**3
  call check_close(d,expected,2e-13_dp,'multinomial density',fails)
  call rmultinomial(100,[1.0_dp,2.0_dp,7.0_dp],counts)
  if (sum(counts)/=100 .or. any(counts<0)) then
    print '(a)','FAIL multinomial sample'; fails=fails+1
  end if

  mu2=0.0_dp; sig2=0.0_dp; sig2(1,1)=1.0_dp; sig2(2,2)=1.0_dp
  d=dmultinormal([1.0_dp,0.0_dp],mu2,sig2)
  expected=exp(-0.5_dp)/(2.0_dp*pi)
  call check_close(d,expected,3e-12_dp,'multinormal density',fails)
  call rmultinormal(8,mu2,sig2,xmv,123)
  if (size(xmv,1)/=8 .or. size(xmv,2)/=2) then
    print '(a)','FAIL multinormal dimensions'; fails=fails+1
  end if

  if(fails/=0) then
    print '(a,i0)','test_distributions: FAIL ',fails
    error stop 1
  end if
  print '(a)','test_distributions: PASS'
contains
  subroutine check_close(a,b,tol,label,fails)
    real(dp),intent(in)::a,b,tol
    character(len=*),intent(in)::label
    integer,intent(inout)::fails
    if(abs(a-b)>tol*max(1.0_dp,abs(a),abs(b)))then
      print '(a,2(1x,es24.16))','FAIL '//trim(label)//':',a,b
      fails=fails+1
    end if
  end subroutine check_close
end program test_distributions
