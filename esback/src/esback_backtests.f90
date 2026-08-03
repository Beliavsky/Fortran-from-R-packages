! SPDX-License-Identifier: GPL-3.0-only
module esback_backtests
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
  use esback_kinds, only: dp
  use esback_types
  use esback_math
  use esback_esreg
  implicit none
  private
  public :: er_backtest, cc_backtest, esr_backtest
contains
  subroutine er_backtest(r,q,e,result,s,b)
    real(dp),intent(in)::r(:),q(:),e(:)
    type(er_backtest_result),intent(out)::result
    real(dp),intent(in),optional::s(:)
    integer,intent(in),optional::b
    real(dp),allocatable::x(:),xs(:)
    integer::nb,n
    real(dp)::p2,p1,nanv
    nanv=ieee_value(1.0_dp,ieee_quiet_nan);result%pvalue_twosided_standardized=nanv;result%pvalue_onesided_standardized=nanv
    if(size(q)/=size(r).or.size(e)/=size(r))then;result%status=esback_invalid_input;return;end if
    nb=1000;if(present(b))nb=b;if(nb<1)then;result%status=esback_invalid_input;return;end if
    n=count(r<=q);result%n_exceedances=n
    if(n<2)then;result%status=esback_insufficient_data;return;end if
    allocate(x(n));x=pack(r-e,r<=q)
    call er_bootstrap_test(x,nb,p2,p1,result%status);if(result%status/=esback_ok)return
    result%pvalue_twosided_simple=p2;result%pvalue_onesided_simple=p1
    if(present(s))then
      if(size(s)/=size(r).or.any(s<=0.0_dp))then;result%status=esback_invalid_input;return;end if
      allocate(xs(n));xs=pack((r-e)/s,r<=q)
      call er_bootstrap_test(xs,nb,p2,p1,result%status);if(result%status/=esback_ok)return
      result%pvalue_twosided_standardized=p2;result%pvalue_onesided_standardized=p1;result%has_standardized=.true.
    end if
  end subroutine er_backtest

  subroutine er_bootstrap_test(x,b,p2,p1,status)
    real(dp),intent(in)::x(:);integer,intent(in)::b
    real(dp),intent(out)::p2,p1;integer,intent(out)::status
    type(rng_state)::rng
    real(dp),allocatable::xb(:),t(:)
    real(dp)::t0,tm,sd
    integer::i,j,m
    sd=sample_sd(x);if(.not.ieee_is_finite(sd).or.sd<=0.0_dp)then;status=esback_insufficient_data;return;end if
    t0=sample_mean(x)/sd*sqrt(real(size(x),dp));allocate(xb(size(x)),t(b));m=0;call rng_seed(rng,1_8)
    do i=1,b
      do j=1,size(x);xb(j)=x(rng_index(rng,size(x)));end do
      sd=sample_sd(xb)
      if(ieee_is_finite(sd).and.sd>0.0_dp)then;m=m+1;t(m)=sample_mean(xb)/sd*sqrt(real(size(x),dp));end if
    end do
    if(m==0)then;status=esback_bootstrap_failed;return;end if
    tm=sum(t(:m))/real(m,dp)
    p2=real(count(abs(t(:m)-tm)>=abs(t0)),dp)/real(m,dp)
    p1=real(count(t(:m)-tm<=t0),dp)/real(m,dp);status=esback_ok
  end subroutine er_bootstrap_test

  subroutine cc_backtest(r,q,e,alpha,result,s,hommel)
    real(dp),intent(in)::r(:),q(:),e(:),alpha
    type(cc_backtest_result),intent(out)::result
    real(dp),intent(in),optional::s(:)
    logical,intent(in),optional::hommel
    real(dp),allocatable::v(:,:),omega(:,:),oinv(:,:),m(:),z(:),p(:),hv(:,:),bvec(:)
    real(dp)::tstat,nanv
    integer::n,st,i
    logical::hm
    nanv=ieee_value(1.0_dp,ieee_quiet_nan);result%pvalue_twosided_general=nanv;result%pvalue_onesided_general=nanv
    n=size(r);hm=.true.;if(present(hommel))hm=hommel
    if(n<3.or.size(q)/=n.or.size(e)/=n.or.alpha<=0.0_dp.or.alpha>=1.0_dp)then;result%status=esback_invalid_input;return;end if
    allocate(v(n,2));v(:,1)=alpha-merge(1.0_dp,0.0_dp,r<=q);v(:,2)=e-q+merge(q-r,0.0_dp,r<=q)/alpha
    omega=matmul(transpose(v),v)/real(n,dp);m=sum(v,dim=1)/real(n,dp);call invert_matrix(omega,oinv,st)
    if(st/=esback_ok)then;result%status=st;return;end if
    tstat=real(n,dp)*dot_product(m,matmul(oinv,m));result%pvalue_twosided_simple=chi_square_survival(tstat,2)
    allocate(z(2),p(2));do i=1,2;z(i)=sqrt(real(n,dp))*m(i)/sqrt(max(omega(i,i),tiny(1.0_dp)));p(i)=1.0_dp-normal_cdf(z(i));end do
    result%pvalue_onesided_simple=multiple_one_sided_pvalue(p,hm)
    if(present(s))then
      if(size(s)/=n.or.any(s<=0.0_dp))then;result%status=esback_invalid_input;return;end if
      allocate(hv(n,4),bvec(n));bvec=((q-e)/(alpha*s))*v(:,1)+(1.0_dp/s)*v(:,2)
      tstat=real(n,dp)*(sample_mean(bvec)**2)/(sum(bvec*bvec)/real(n,dp))
      result%pvalue_twosided_general=chi_square_survival(tstat,1)
      hv(:,1)=v(:,1);hv(:,2)=abs(q)*v(:,1);hv(:,3)=v(:,2);hv(:,4)=v(:,2)/s
      deallocate(z,p);allocate(z(4),p(4))
      do i=1,4
        z(i)=sqrt(real(n,dp))*sample_mean(hv(:,i))/sqrt(max(sum(hv(:,i)**2)/real(n,dp),tiny(1.0_dp)))
        p(i)=1.0_dp-normal_cdf(z(i))
      end do
      result%pvalue_onesided_general=multiple_one_sided_pvalue(p,hm);result%has_general=.true.
    end if
    result%status=esback_ok
  end subroutine cc_backtest

  real(dp) function multiple_one_sided_pvalue(p,hommel) result(pv)
    real(dp),intent(in)::p(:);logical,intent(in)::hommel
    real(dp),allocatable::ps(:);real(dp)::hm
    integer::i,m
    m=size(p);allocate(ps(m));ps=p;call sort_real(ps)
    if(hommel)then
      hm=sum([(1.0_dp/real(i,dp),i=1,m)])
      pv=min(1.0_dp,real(m,dp)*hm*minval(ps/[(real(i,dp),i=1,m)]))
    else
      pv=min(1.0_dp,real(m,dp)*minval(ps))
    end if
  end function multiple_one_sided_pvalue

  subroutine esr_backtest(r,e,alpha,version,result,q,b,options)
    real(dp),intent(in)::r(:),e(:),alpha
    integer,intent(in)::version
    type(esr_backtest_result),intent(out)::result
    real(dp),intent(in),optional::q(:)
    integer,intent(in),optional::b
    type(esreg_options),intent(in),optional::options
    type(esreg_options)::opt
    real(dp),allocatable::xq(:,:),xe(:,:),y(:),svec(:),subcov(:,:),invsub(:,:),sb(:),tb(:)
    integer,allocatable::idx(:)
    type(esreg_fit_result)::fitb
    type(rng_state)::rng
    integer::n,nb,kq,ke,st,i,j,success,fail
    real(dp)::t0,se,nanv
    opt=esreg_options();if(present(options))opt=options;nb=0;if(present(b))nb=b;n=size(r);nanv=ieee_value(1.0_dp,ieee_quiet_nan)
    result%version = version
    result%pvalue_onesided_asymptotic = nanv
    result%pvalue_onesided_bootstrap = nanv
    result%pvalue_twosided_bootstrap = nanv
    if (size(e) /= n .or. n < 10 .or. alpha <= 0.0_dp .or. &
        alpha >= 1.0_dp .or. version < 1 .or. version > 3) then
      result%status = esback_invalid_input
      return
    end if
    select case(version)
    case(1)
      kq=2;ke=2;allocate(xq(n,kq),xe(n,ke),y(n));xq(:,1)=1.0_dp;xq(:,2)=e;xe=xq;y=r
    case(2)
      if(.not.present(q))then;result%status=esback_invalid_input;return;end if
      if(size(q)/=n)then;result%status=esback_invalid_input;return;end if
      kq=2;ke=2;allocate(xq(n,kq),xe(n,ke),y(n));xq(:,1)=1.0_dp;xq(:,2)=q;xe(:,1)=1.0_dp;xe(:,2)=e;y=r
    case(3)
      kq=2;ke=1;allocate(xq(n,kq),xe(n,ke),y(n));xq(:,1)=1.0_dp;xq(:,2)=e;xe(:,1)=1.0_dp;y=r-e;result%one_sided=.true.
    end select
    call esreg_fit(xq, xe, y, alpha, result%fit, opt, .true.)
    if (result%fit%status /= esback_ok .or. &
        .not. result%fit%covariance_available) then
      result%status = result%fit%status
      return
    end if
    if(version<=2)then
      allocate(svec(2),subcov(2,2));svec=[result%fit%coefficients_e(1),result%fit%coefficients_e(2)-1.0_dp]
      subcov = result%fit%covariance(kq+1:kq+2,kq+1:kq+2)
      call invert_matrix(subcov, invsub, st)
      if (st /= esback_ok) then
        result%status = st
        return
      end if
      t0=dot_product(svec,matmul(invsub,svec));result%statistic=t0;result%pvalue_twosided_asymptotic=chi_square_survival(t0,2)
    else
      se=sqrt(max(result%fit%covariance(kq+1,kq+1),tiny(1.0_dp)));t0=result%fit%coefficients_e(1)/se;result%statistic=t0
      result%pvalue_onesided_asymptotic=normal_cdf(t0);result%pvalue_twosided_asymptotic=2.0_dp*(1.0_dp-normal_cdf(abs(t0)))
    end if
    if(nb>0)then
      allocate(idx(n),tb(nb));success=0;fail=0;call rng_seed(rng,opt%seed+7919_8)
      do i=1,nb
        do j=1,n;idx(j)=rng_index(rng,n);end do
        call esreg_fit(xq(idx,:),xe(idx,:),y(idx),alpha,fitb,opt,.true.)
        if(fitb%status/=esback_ok.or..not.fitb%covariance_available)then;fail=fail+1;cycle;end if
        if(version<=2)then
          sb = fitb%coefficients_e - result%fit%coefficients_e
          subcov = fitb%covariance(kq+1:kq+2,kq+1:kq+2)
          call invert_matrix(subcov, invsub, st)
          if(st/=esback_ok)then;fail=fail+1;cycle;end if
          success=success+1;tb(success)=dot_product(sb,matmul(invsub,sb))
        else
          se = sqrt(max(fitb%covariance(kq+1,kq+1), tiny(1.0_dp)))
          success = success + 1
          tb(success) = (fitb%coefficients_e(1) - &
            result%fit%coefficients_e(1))/se
        end if
      end do
      if(real(fail,dp)/real(nb,dp)>=0.05_dp.or.success==0)then;result%status=esback_bootstrap_failed;return;end if
      result%bootstrap_used=.true.;result%bootstrap_successes=success
      if(version<=2)then
        result%pvalue_twosided_bootstrap=real(count(tb(:success)>=t0),dp)/real(success,dp)
      else
        result%pvalue_twosided_bootstrap=real(count(abs(tb(:success))>=abs(t0)),dp)/real(success,dp)
        result%pvalue_onesided_bootstrap=real(count(tb(:success)<=t0),dp)/real(success,dp)
      end if
    end if
    result%status=esback_ok
  end subroutine esr_backtest
end module esback_backtests
