program test_sampling_dependence
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use mc2d, only : dp, mcnode, mcstoc, lhs, rtrunc, cornode, correlation_dp
  implicit none
  integer,parameter::n=1000
  real(dp)::data(n,2),target(2,2),rho
  real(dp),allocatable::correlated(:,:)
  real(dp)::lh(10,5,2),tr(500)
  type(mcnode)::xv
  integer::i,j,k,fails

  fails=0
  xv=mcstoc(uniform_sampler,type='V',nsv=25,seed=11)
  if(any(shape(xv%value)/=[25,1,1]))call fail('mcstoc dimensions',fails)
  if(any(xv%value<0.0_dp).or.any(xv%value>1.0_dp))call fail('mcstoc values',fails)

  call lhs(identity_quantile,10,5,2,lh)
  if(any(lh<0.0_dp).or.any(lh>1.0_dp))call fail('lhs support',fails)
  do k=1,2;do j=1,5
    if(.not.strata_once(lh(:,j,k)))call fail('lhs strata',fails)
  end do;end do

  call rtrunc(normal_cdf_local,normal_quantile_local,500,-1.0_dp,2.0_dp,tr)
  if(any(ieee_is_nan(tr)))call fail('rtrunc NaN',fails)
  if(any(tr<=-1.0_dp).or.any(tr>2.0_dp))call fail('rtrunc support',fails)

  do i=1,n
    data(i,1)=real(i,dp)
    data(i,2)=real(mod(7919*i,1009),dp)+0.001_dp*real(i,dp)
  end do
  target=reshape([1.0_dp,0.8_dp,0.8_dp,1.0_dp],[2,2])
  correlated=cornode(data,target,seed=42)
  if(any(correlated(:,1)/=data(:,1)))call fail('cornode preserves first column',fails)
  rho=correlation_dp(correlated(:,1),correlated(:,2),'spearman')
  if(abs(rho-0.8_dp)>0.05_dp)then
    print '(a,f10.6)','FAIL cornode rho = ',rho;fails=fails+1
  end if

  if(fails/=0)then
    print '(a,i0)','test_sampling_dependence: FAIL ',fails;error stop 1
  end if
  print '(a,f10.6)','test_sampling_dependence: PASS, rho=',rho
contains
  subroutine uniform_sampler(n,x)
    integer,intent(in)::n;real(dp),intent(out)::x(n);call random_number(x)
  end subroutine
  real(dp) function identity_quantile(p) result(q)
    real(dp),intent(in)::p;q=p
  end function
  real(dp) function normal_cdf_local(x) result(p)
    use mvtnorm_special,only:normal_cdf
    real(dp),intent(in)::x;p=normal_cdf(x)
  end function
  real(dp) function normal_quantile_local(p) result(q)
    use mvtnorm_special,only:normal_quantile
    real(dp),intent(in)::p;q=normal_quantile(p)
  end function
  logical function strata_once(x) result(ok)
    real(dp),intent(in)::x(:);integer::bins(size(x)),ii
    bins=1+int(x*real(size(x),dp));bins=min(bins,size(x));ok=.true.
    do ii=1,size(x)
      if(count(bins==ii)/=1)then;ok=.false.;return;end if
    end do
  end function
  subroutine fail(label,fails)
    character(len=*),intent(in)::label;integer,intent(inout)::fails
    print '(a)','FAIL '//trim(label);fails=fails+1
  end subroutine
end program test_sampling_dependence
