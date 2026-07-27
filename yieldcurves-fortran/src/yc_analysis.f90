! SPDX-License-Identifier: MIT
! Copyright (c) 2026 Charles Coverdale
module yc_analysis
  use yc_kinds, only : dp
  use yc_types, only : curve_t, series_t, duration_result_t, bond_duration_result_t
  use yc_types, only : zspread_result_t, carry_result_t, slope_result_t, factor_result_t
  use yc_utils, only : valid_positive_vector, valid_finite_vector, sort_pairs, linear_interpolate
  use yc_utils, only : lower_string, quiet_nan, payment_times
  use yc_curve_ops, only : yc_predict
  implicit none
  private
  public :: yc_par_to_zero, yc_zero_to_par, yc_duration, yc_bond_duration
  public :: yc_zspread, yc_key_rate_duration, yc_carry, yc_slope
  public :: yc_level_slope_curvature

contains

  function yc_par_to_zero(maturities,par_rates,frequency) result(out)
    real(dp),intent(in)::maturities(:),par_rates(:)
    integer,intent(in),optional::frequency
    type(series_t)::out
    real(dp),allocatable::m(:),p(:),coupon_dates(:),intermediate(:),zint(:)
    real(dp)::coupon,pv_coupons,pv_final,periods
    integer::freq,i,j,ncf
    logical::ok
    freq=1;if(present(frequency))freq=frequency
    if(size(maturities)/=size(par_rates).or..not.valid_positive_vector(maturities).or. &
      .not.valid_finite_vector(par_rates).or.(freq/=1.and.freq/=2))then
      out%ok=.false.;out%message='Invalid par-to-zero inputs.';return
    end if
    allocate(m,source=maturities);allocate(p,source=par_rates);call sort_pairs(m,p)
    allocate(out%x,source=m);allocate(out%y(size(m)));out%y=0.0_dp
    out%y(1)=p(1)
    do i=2,size(m)
      coupon_dates=payment_times(m(i),freq,ok)
      if(.not.ok)then;out%ok=.false.;out%message='Maturity must align with coupon frequency.';return;end if
      ncf=size(coupon_dates)
      if(ncf<=1)then;out%y(i)=p(i);cycle;end if
      coupon=p(i)/real(freq,dp)
      allocate(intermediate(ncf-1),zint(ncf-1))
      intermediate=coupon_dates(1:ncf-1)
      do j=1,ncf-1
        zint(j)=linear_interpolate(m(1:i-1),out%y(1:i-1),intermediate(j))
      end do
      pv_coupons=sum(coupon/(1.0_dp+zint/real(freq,dp))**(intermediate*real(freq,dp)))
      pv_final=1.0_dp-pv_coupons
      periods=m(i)*real(freq,dp)
      if(pv_final<=0.0_dp)then;out%ok=.false.;out%message='Bootstrap produced non-positive final PV.';return;end if
      out%y(i)=real(freq,dp)*(((1.0_dp+coupon)/pv_final)**(1.0_dp/periods)-1.0_dp)
      deallocate(intermediate,zint)
    end do
  end function yc_par_to_zero

  function yc_zero_to_par(maturities,zero_rates,frequency) result(out)
    real(dp),intent(in)::maturities(:),zero_rates(:)
    integer,intent(in),optional::frequency
    type(series_t)::out
    real(dp),allocatable::m(:),z(:),times(:),zat(:),df(:)
    integer::freq,i,j
    logical::ok
    freq=1;if(present(frequency))freq=frequency
    if(size(maturities)/=size(zero_rates).or..not.valid_positive_vector(maturities).or. &
      .not.valid_finite_vector(zero_rates).or.(freq/=1.and.freq/=2))then
      out%ok=.false.;out%message='Invalid zero-to-par inputs.';return
    end if
    allocate(m,source=maturities);allocate(z,source=zero_rates);call sort_pairs(m,z)
    allocate(out%x,source=m);allocate(out%y(size(m)))
    do i=1,size(m)
      times=payment_times(m(i),freq,ok)
      if(.not.ok)then;out%ok=.false.;out%message='Maturity must align with coupon frequency.';return;end if
      if(size(times)<=1)then;out%y(i)=z(i);cycle;end if
      allocate(zat(size(times)),df(size(times)))
      do j=1,size(times);zat(j)=linear_interpolate(m(1:i),z(1:i),times(j));end do
      df=(1.0_dp+zat/real(freq,dp))**(-times*real(freq,dp))
      out%y(i)=real(freq,dp)*(1.0_dp-df(size(df)))/sum(df)
      deallocate(zat,df)
    end do
  end function yc_zero_to_par

  function yc_duration(curve,maturities,compounding) result(out)
    type(curve_t),intent(in)::curve
    real(dp),intent(in),optional::maturities(:)
    character(len=*),intent(in),optional::compounding
    type(duration_result_t)::out
    type(series_t)::pred
    character(len=20)::comp
    real(dp),allocatable::m(:)
    comp='continuous';if(present(compounding))comp=trim(lower_string(compounding))
    if(present(maturities))then;allocate(m,source=maturities);else;allocate(m,source=curve%maturities);end if
    pred=yc_predict(curve,m)
    if(.not.pred%ok)then;out%ok=.false.;out%message=pred%message;return;end if
    allocate(out%maturity,source=m);allocate(out%macaulay(size(m)),out%modified(size(m)),out%convexity(size(m)))
    out%macaulay=m
    select case(comp)
    case('continuous')
      out%modified=m;out%convexity=m*m
    case('annual')
      out%modified=m/(1.0_dp+pred%y)
      out%convexity=m*(m+1.0_dp)/(1.0_dp+pred%y)**2
    case('semi_annual')
      out%modified=m/(1.0_dp+pred%y/2.0_dp)
      out%convexity=m*(m+0.5_dp)/(1.0_dp+pred%y/2.0_dp)**2
    case default
      out%ok=.false.;out%message='Unknown compounding convention.'
    end select
  end function yc_duration

  function yc_bond_duration(face,coupon_rate,maturity,yield_rate,frequency,compounding) result(out)
    real(dp),intent(in)::face,coupon_rate,maturity,yield_rate
    integer,intent(in),optional::frequency
    character(len=*),intent(in),optional::compounding
    type(bond_duration_result_t)::out
    real(dp),allocatable::times(:),cfs(:),disc(:)
    real(dp)::yper
    integer::freq,n
    logical::ok
    character(len=20)::comp
    freq=2;if(present(frequency))freq=frequency
    comp='semi_annual';if(present(compounding))comp=trim(lower_string(compounding))
    if(face<=0.0_dp.or.maturity<=0.0_dp.or.(freq/=1.and.freq/=2))then
      out%ok=.false.;out%message='Invalid bond-duration inputs.';return
    end if
    if(comp=='annual')freq=1
    times=payment_times(maturity,freq,ok)
    if(.not.ok)then;out%ok=.false.;out%message='Maturity must align with coupon frequency.';return;end if
    n=size(times);allocate(cfs(n),disc(n));cfs=face*coupon_rate/real(freq,dp);cfs(n)=cfs(n)+face
    if(comp=='continuous')then
      disc=exp(-yield_rate*times)
      out%price=sum(cfs*disc)
      out%macaulay_duration=sum(times*cfs*disc)/out%price
      out%modified_duration=out%macaulay_duration
      out%convexity=sum(times*times*cfs*disc)/out%price
    else if(comp=='annual'.or.comp=='semi_annual')then
      yper=yield_rate/real(freq,dp)
      if(1.0_dp+yper<=0.0_dp)then;out%ok=.false.;out%message='Invalid yield.';return;end if
      disc=(1.0_dp+yper)**(-times*real(freq,dp))
      out%price=sum(cfs*disc)
      out%macaulay_duration=sum(times*cfs*disc)/out%price
      out%modified_duration=out%macaulay_duration/(1.0_dp+yper)
      out%convexity=sum(times*(times+1.0_dp/real(freq,dp))*cfs*disc)/(out%price*(1.0_dp+yper)**2)
    else
      out%ok=.false.;out%message='Unknown compounding convention.'
    end if
  end function yc_bond_duration

  function yc_zspread(price,coupon_rate,maturity,curve,face,frequency) result(out)
    real(dp),intent(in)::price,coupon_rate,maturity
    type(curve_t),intent(in)::curve
    real(dp),intent(in),optional::face
    integer,intent(in),optional::frequency
    type(zspread_result_t)::out
    real(dp)::a,b,fa,fb,mid,fm,facev
    real(dp),allocatable::times(:),cfs(:),z(:)
    integer::freq,i,iter
    logical::ok
    facev=100.0_dp;if(present(face))facev=face
    freq=2;if(present(frequency))freq=frequency
    out%price=price
    if(price<=0.0_dp.or.facev<=0.0_dp.or.maturity<=0.0_dp.or.(freq/=1.and.freq/=2))then
      out%ok=.false.;out%message='Invalid z-spread inputs.';return
    end if
    times=payment_times(maturity,freq,ok)
    if(.not.ok)then;out%ok=.false.;out%message='Maturity must align with coupon frequency.';return;end if
    allocate(cfs(size(times)),z(size(times)))
    cfs=facev*coupon_rate/real(freq,dp);cfs(size(cfs))=cfs(size(cfs))+facev
    do i=1,size(times);z(i)=benchmark_rate_at(curve,times(i));end do
    a=-0.05_dp;b=0.50_dp;fa=bond_error(a);fb=bond_error(b)
    do iter=1,60
      if(fa*fb<=0.0_dp)exit
      a=max(-0.99_dp-minval(z)+1.0e-10_dp,a-0.05_dp)
      b=b*1.5_dp+0.05_dp
      fa=bond_error(a);fb=bond_error(b)
    end do
    if(fa*fb>0.0_dp)then;out%ok=.false.;out%message='Could not bracket z-spread.';return;end if
    do iter=1,200
      mid=0.5_dp*(a+b);fm=bond_error(mid)
      if(abs(fm)<1.0e-12_dp.or.abs(b-a)<1.0e-12_dp)exit
      if(fa*fm<=0.0_dp)then;b=mid;fb=fm;else;a=mid;fa=fm;end if
    end do
    out%zspread=mid;out%model_price=price+fm;out%iterations=iter
  contains
    real(dp) function bond_error(s) result(v)
      real(dp),intent(in)::s
      if(any(1.0_dp+z+s<=0.0_dp))then;v=huge(1.0_dp);return;end if
      v=sum(cfs*(1.0_dp+z+s)**(-times))-price
    end function bond_error
  end function yc_zspread

  function yc_key_rate_duration(coupon_rate,maturity,curve,key_rates,shift,face,frequency) result(out)
    real(dp),intent(in)::coupon_rate,maturity
    type(curve_t),intent(in)::curve
    real(dp),intent(in),optional::key_rates(:),shift,face
    integer,intent(in),optional::frequency
    type(series_t)::out
    real(dp),allocatable::kr(:),times(:),cfs(:),zbase(:),zbump(:),bump(:)
    real(dp)::shiftv,facev,base_price,price_up,left,right,t
    integer::freq,i,j
    logical::ok
    real(dp),parameter::default_kr(5)=[1.0_dp,2.0_dp,5.0_dp,10.0_dp,30.0_dp]
    shiftv=0.0001_dp;if(present(shift))shiftv=shift
    facev=100.0_dp;if(present(face))facev=face
    freq=2;if(present(frequency))freq=frequency
    if(present(key_rates))then;allocate(kr,source=key_rates);else;allocate(kr,source=default_kr);end if
    if(maturity<=0.0_dp.or.shiftv<=0.0_dp.or.facev<=0.0_dp.or. &
      .not.valid_positive_vector(kr).or.(freq/=1.and.freq/=2))then
      out%ok=.false.;out%message='Invalid key-rate-duration inputs.';return
    end if
    call sort_one(kr)
    times=payment_times(maturity,freq,ok)
    if(.not.ok)then;out%ok=.false.;out%message='Maturity must align with coupon frequency.';return;end if
    allocate(cfs(size(times)),zbase(size(times)),zbump(size(times)),bump(size(times)))
    cfs=facev*coupon_rate/real(freq,dp);cfs(size(cfs))=cfs(size(cfs))+facev
    do j=1,size(times);zbase(j)=benchmark_rate_at(curve,times(j));end do
    base_price=sum(cfs*(1.0_dp+zbase)**(-times))
    allocate(out%x,source=kr);allocate(out%y(size(kr)))
    do i=1,size(kr)
      if(i>1)then;left=kr(i-1);else;left=0.0_dp;end if
      if(i<size(kr))then;right=kr(i+1);else;right=kr(i)+(kr(i)-left);end if
      bump=0.0_dp
      do j=1,size(times)
        t=times(j)
        if(t>=left.and.t<=kr(i))then
          if(kr(i)>left)bump(j)=(t-left)/(kr(i)-left)
        else if(t>kr(i).and.t<=right)then
          if(right>kr(i))bump(j)=(right-t)/(right-kr(i))
        end if
      end do
      zbump=zbase+bump*shiftv
      price_up=sum(cfs*(1.0_dp+zbump)**(-times))
      out%y(i)=-(price_up-base_price)/(shiftv*base_price)
    end do
  end function yc_key_rate_duration

  function yc_carry(curve,maturities,horizon,funding_rate) result(out)
    type(curve_t),intent(in)::curve
    real(dp),intent(in),optional::maturities(:),horizon,funding_rate
    type(carry_result_t)::out
    real(dp),allocatable::m(:)
    real(dp)::h,fr
    type(series_t)::current,rolled,fund
    integer::i,n
    h=1.0_dp/12.0_dp;if(present(horizon))h=horizon
    if(h<=0.0_dp)then;out%ok=.false.;out%message='Horizon must be positive.';return;end if
    if(present(maturities))then
      allocate(m,source=maturities)
    else
      n=count(curve%maturities>h);allocate(m(n));n=0
      do i=1,size(curve%maturities);if(curve%maturities(i)>h)then;n=n+1;m(n)=curve%maturities(i);end if;end do
    end if
    if(size(m)==0.or.any(m<=h))then;out%ok=.false.;out%message='Maturities must exceed horizon.';return;end if
    if(present(funding_rate))then;fr=funding_rate;else;fund=yc_predict(curve,[h]);fr=fund%y(1);end if
    current=yc_predict(curve,m);rolled=yc_predict(curve,m-h)
    allocate(out%maturity,source=m);allocate(out%carry(size(m)),out%rolldown(size(m)),out%total(size(m)))
    out%carry=(current%y-fr)*h
    out%rolldown=(current%y-rolled%y)*(m-h)
    out%total=out%carry+out%rolldown
  end function yc_carry

  function yc_slope(curve) result(out)
    type(curve_t),intent(in)::curve
    type(slope_result_t)::out
    real(dp)::r3m,r2,r5,r10,r30
    r3m=safe_rate(curve,0.25_dp);r2=safe_rate(curve,2.0_dp);r5=safe_rate(curve,5.0_dp)
    r10=safe_rate(curve,10.0_dp);r30=safe_rate(curve,30.0_dp)
    out%spread_2s10s=r10-r2;out%spread_2s30s=r30-r2;out%spread_5s30s=r30-r5
    out%spread_3m10y=r10-r3m;out%butterfly_2s5s10s=2.0_dp*r5-r2-r10
  end function yc_slope

  function yc_level_slope_curvature(curve) result(out)
    type(curve_t),intent(in)::curve
    type(factor_result_t)::out
    real(dp),allocatable::r(:)
    integer::mid
    if(trim(curve%method)=='nelson_siegel'.or.trim(curve%method)=='svensson')then
      out%level=curve%beta0;out%slope=curve%beta1;out%curvature=curve%beta2
    else
      if(allocated(curve%fitted))then;allocate(r,source=curve%fitted);else;allocate(r,source=curve%rates);end if
      mid=(size(r)+1)/2
      out%level=sum(r)/real(size(r),dp);out%slope=r(1)-r(size(r));out%curvature=2.0_dp*r(mid)-r(1)-r(size(r))
    end if
  end function yc_level_slope_curvature

  real(dp) function benchmark_rate_at(curve,m) result(r)
    type(curve_t),intent(in)::curve
    real(dp),intent(in)::m
    if(allocated(curve%fitted))then
      r=linear_interpolate(curve%maturities,curve%fitted,m)
    else
      r=linear_interpolate(curve%maturities,curve%rates,m)
    end if
  end function benchmark_rate_at

  real(dp) function curve_rate_at(curve,m) result(r)
    type(curve_t),intent(in)::curve
    real(dp),intent(in)::m
    type(series_t)::p
    p=yc_predict(curve,[m]);r=p%y(1)
  end function curve_rate_at

  real(dp) function safe_rate(curve,m) result(r)
    type(curve_t),intent(in)::curve
    real(dp),intent(in)::m
    if(m<minval(curve%maturities).or.m>maxval(curve%maturities))then;r=quiet_nan();else;r=curve_rate_at(curve,m);end if
  end function safe_rate

  subroutine sort_one(x)
    real(dp),intent(inout)::x(:)
    real(dp)::v
    integer::i,j
    do i=2,size(x);v=x(i);j=i-1;do while(j>=1);if(x(j)<=v)exit;x(j+1)=x(j);j=j-1;end do;x(j+1)=v;end do
  end subroutine sort_one

end module yc_analysis
