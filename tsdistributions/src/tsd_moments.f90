! SPDX-License-Identifier: GPL-2.0-only
module tsd_moments
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
  use ghyp_kinds, only : dp, pi
  use ghyp_gig, only : gig_raw_moment
  use tsd_types, only : moment_summary, authorized_domain_result, tsd_success, &
                        tsd_invalid_argument, distribution_id, dist_norm, dist_std, &
                        dist_snorm, dist_sstd, dist_ged, dist_sged, dist_nig, &
                        dist_gh, dist_jsu, dist_ghst
  use tsd_distributions, only : paramgh, paramghst
  implicit none
  private

  public :: dskewness, dkurtosis, distribution_moments, authorized_domain

contains

  real(dp) function dskewness(distribution, skew, shape, lambda) result(value)
    character(len=*), intent(in) :: distribution
    real(dp), intent(in), optional :: skew, shape, lambda
    real(dp) :: xi, nu, lam
    xi=1.0_dp;nu=5.0_dp;lam=-0.5_dp
    if(present(skew))xi=skew;if(present(shape))nu=shape;if(present(lambda))lam=lambda
    select case(distribution_id(distribution))
    case(dist_norm,dist_std,dist_ged)
      value=0.0_dp
    case(dist_snorm)
      value=snorm_skewness(xi)
    case(dist_sstd)
      value=sstd_skewness(xi,nu)
    case(dist_sged)
      value=sged_skewness(xi,nu)
    case(dist_nig)
      value=gh_skewness(xi,nu,-0.5_dp)
    case(dist_gh)
      value=gh_skewness(xi,nu,lam)
    case(dist_jsu)
      value=jsu_skewness(xi,nu)
    case(dist_ghst)
      value=ghst_skewness(xi,nu)
    case default
      value=ieee_value(0.0_dp,ieee_quiet_nan)
    end select
  end function dskewness

  real(dp) function dkurtosis(distribution, skew, shape, lambda) result(value)
    character(len=*), intent(in) :: distribution
    real(dp), intent(in), optional :: skew, shape, lambda
    real(dp) :: xi, nu, lam
    xi=1.0_dp;nu=5.0_dp;lam=-0.5_dp
    if(present(skew))xi=skew;if(present(shape))nu=shape;if(present(lambda))lam=lambda
    select case(distribution_id(distribution))
    case(dist_norm,dist_snorm)
      value=0.0_dp
    case(dist_std)
      if(nu>4.0_dp)then;value=6.0_dp/(nu-4.0_dp);else;value=ieee_value(0.0_dp,ieee_quiet_nan);end if
    case(dist_sstd)
      value=sstd_excess_kurtosis(xi,nu)
    case(dist_ged)
      value=ged_excess_kurtosis(nu)
    case(dist_sged)
      value=sged_excess_kurtosis(xi,nu)
    case(dist_nig)
      value=gh_excess_kurtosis(xi,nu,-0.5_dp)
    case(dist_gh)
      value=gh_excess_kurtosis(xi,nu,lam)
    case(dist_jsu)
      value=jsu_excess_kurtosis(xi,nu)
    case(dist_ghst)
      value=ghst_excess_kurtosis(xi,nu)
    case default
      value=ieee_value(0.0_dp,ieee_quiet_nan)
    end select
  end function dkurtosis

  function distribution_moments(distribution, mu, sigma, skew, shape, lambda) result(summary)
    character(len=*),intent(in)::distribution
    real(dp),intent(in),optional::mu,sigma,skew,shape,lambda
    type(moment_summary)::summary
    real(dp)::m,s,xi,nu,lam
    m=0.0_dp;s=1.0_dp;xi=1.0_dp;nu=5.0_dp;lam=-0.5_dp
    if(present(mu))m=mu;if(present(sigma))s=sigma;if(present(skew))xi=skew;if(present(shape))nu=shape;if(present(lambda))lam=lambda
    if(s<=0.0_dp)then;summary%status=tsd_invalid_argument;summary%message='sigma must be positive';return;end if
    summary%mean=m;summary%standard_deviation=s
    summary%skewness=dskewness(distribution,xi,nu,lam)
    summary%excess_kurtosis=dkurtosis(distribution,xi,nu,lam)
    if(.not.ieee_is_finite(summary%skewness).or..not.ieee_is_finite(summary%excess_kurtosis))then
      summary%status=tsd_invalid_argument;summary%message='requested moment does not exist'
    else
      summary%status=tsd_success;summary%message=''
    end if
  end function distribution_moments

  pure real(dp) function snorm_skewness(xi) result(value)
    real(dp),intent(in)::xi
    real(dp)::m1,m2,m3
    m1=sqrt(2.0_dp/pi);m2=1.0_dp;m3=2.0_dp*sqrt(2.0_dp/pi)
    value=(xi-1.0_dp/xi)*((m3+2.0_dp*m1**3-3.0_dp*m1*m2)*(xi*xi+xi**(-2))+3.0_dp*m1*m2-4.0_dp*m1**3)/ &
      (((m2-m1*m1)*(xi*xi+xi**(-2))+2.0_dp*m1*m1-m2)**1.5_dp)
  end function snorm_skewness

  pure real(dp) function sstd_skewness(xi,eta) result(value)
    real(dp),intent(in)::xi,eta
    real(dp)::lambda_s,cx,a,b,my2,my3
    if(eta<=3.0_dp)then;value=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    lambda_s=(xi*xi-1.0_dp)/(xi*xi+1.0_dp)
    cx=exp(log_gamma(0.5_dp*(eta+1.0_dp))-log_gamma(0.5_dp*eta)-0.5_dp*log(pi*(eta-2.0_dp)))
    a=4.0_dp*lambda_s*cx*(eta-2.0_dp)/(eta-1.0_dp);b=sqrt(1.0_dp+3.0_dp*lambda_s**2-a*a)
    my2=1.0_dp+3.0_dp*lambda_s**2
    my3=16.0_dp*cx*lambda_s*(1.0_dp+lambda_s**2)*(eta-2.0_dp)**2/((eta-1.0_dp)*(eta-3.0_dp))
    value=(my3-3.0_dp*a*my2+2.0_dp*a**3)/b**3
  end function sstd_skewness

  pure real(dp) function sstd_excess_kurtosis(xi,eta) result(value)
    real(dp),intent(in)::xi,eta
    real(dp)::lambda_s,cx,a,b,my2,my3,my4
    if(eta<=4.0_dp)then;value=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    lambda_s=(xi*xi-1.0_dp)/(xi*xi+1.0_dp)
    cx=exp(log_gamma(0.5_dp*(eta+1.0_dp))-log_gamma(0.5_dp*eta)-0.5_dp*log(pi*(eta-2.0_dp)))
    a=4.0_dp*lambda_s*cx*(eta-2.0_dp)/(eta-1.0_dp);b=sqrt(1.0_dp+3.0_dp*lambda_s**2-a*a)
    my2=1.0_dp+3.0_dp*lambda_s**2
    my3=16.0_dp*cx*lambda_s*(1.0_dp+lambda_s**2)*(eta-2.0_dp)**2/((eta-1.0_dp)*(eta-3.0_dp))
    my4=3.0_dp*(eta-2.0_dp)*(1.0_dp+10.0_dp*lambda_s**2+5.0_dp*lambda_s**4)/(eta-4.0_dp)
    value=-3.0_dp+(my4-4.0_dp*a*my3+6.0_dp*a*a*my2-3.0_dp*a**4)/b**4
  end function sstd_excess_kurtosis

  pure real(dp) function ged_excess_kurtosis(nu) result(value)
    real(dp),intent(in)::nu
    if(nu<=0.0_dp)then;value=ieee_value(0.0_dp,ieee_quiet_nan);else
      value=exp(2.0_dp*(log_gamma(1.0_dp/nu)-log_gamma(3.0_dp/nu))+log_gamma(5.0_dp/nu)-log_gamma(1.0_dp/nu))-3.0_dp
    end if
  end function ged_excess_kurtosis

  pure real(dp) function sged_skewness(xi,nu) result(value)
    real(dp),intent(in)::xi,nu
    real(dp)::lambda_g,m1,m2,m3
    lambda_g=sqrt(0.5_dp**(2.0_dp/nu)*gamma(1.0_dp/nu)/gamma(3.0_dp/nu))
    m1=2.0_dp**(1.0_dp/nu)*lambda_g*gamma(2.0_dp/nu)/gamma(1.0_dp/nu)
    m2=1.0_dp;m3=(2.0_dp**(1.0_dp/nu)*lambda_g)**3*gamma(4.0_dp/nu)/gamma(1.0_dp/nu)
    value=(xi-1.0_dp/xi)*((m3+2.0_dp*m1**3-3.0_dp*m1*m2)*(xi*xi+xi**(-2))+3.0_dp*m1*m2-4.0_dp*m1**3)/ &
      (((m2-m1*m1)*(xi*xi+xi**(-2))+2.0_dp*m1*m1-m2)**1.5_dp)
  end function sged_skewness

  pure real(dp) function sged_excess_kurtosis(xi,nu) result(value)
    real(dp),intent(in)::xi,nu
    real(dp)::lambda_g,m1,m2,m3,m4,cm4,var2
    lambda_g=sqrt(0.5_dp**(2.0_dp/nu)*gamma(1.0_dp/nu)/gamma(3.0_dp/nu))
    m1=2.0_dp**(1.0_dp/nu)*lambda_g*gamma(2.0_dp/nu)/gamma(1.0_dp/nu);m2=1.0_dp
    m3=(2.0_dp**(1.0_dp/nu)*lambda_g)**3*gamma(4.0_dp/nu)/gamma(1.0_dp/nu)
    m4=(2.0_dp**(1.0_dp/nu)*lambda_g)**4*gamma(5.0_dp/nu)/gamma(1.0_dp/nu)
    cm4=-3.0_dp*m1**4*(xi-1.0_dp/xi)**4+6.0_dp*m1*m1*(xi-1.0_dp/xi)**2*m2*(xi**3+xi**(-3))/(xi+1.0_dp/xi)- &
      4.0_dp*m1*(xi-1.0_dp/xi)*m3*(xi**4-xi**(-4))/(xi+1.0_dp/xi)+m4*(xi**5+xi**(-5))/(xi+1.0_dp/xi)
    var2=((m2-m1*m1)*(xi*xi+xi**(-2))+2.0_dp*m1*m1-m2)**2
    value=cm4/var2-3.0_dp
  end function sged_excess_kurtosis

  pure real(dp) function jsu_skewness(skew,shape) result(value)
    real(dp),intent(in)::skew,shape
    real(dp)::omega,w,s3,var
    omega=-skew/shape;w=exp(shape**(-2));s3=-0.25_dp*sqrt(w)*(w-1.0_dp)**2*(w*(w+2.0_dp)*sinh(3.0_dp*omega)+3.0_dp*sinh(omega))
    var=0.5_dp*(w-1.0_dp)*(w*cosh(2.0_dp*omega)+1.0_dp);value=s3/var**1.5_dp
  end function jsu_skewness

  pure real(dp) function jsu_excess_kurtosis(skew,shape) result(value)
    real(dp),intent(in)::skew,shape
    real(dp)::omega,w,s4,var
    omega=-skew/shape;w=exp(shape**(-2));s4=0.125_dp*(w-1.0_dp)**2*(w*w*(w**4+2.0_dp*w**3+3.0_dp*w*w-3.0_dp)*cosh(4.0_dp*omega)+ &
      4.0_dp*w*w*(w+2.0_dp)*cosh(2.0_dp*omega)+3.0_dp*(2.0_dp*w+1.0_dp))
    var=0.5_dp*(w-1.0_dp)*(w*cosh(2.0_dp*omega)+1.0_dp);value=s4/(var*var)-3.0_dp
  end function jsu_excess_kurtosis

  real(dp) function gh_skewness(skew,shape,lambda) result(value)
    real(dp),intent(in)::skew,shape,lambda
    real(dp)::m2,m3,m4,var,c3,c4
    call gh_raw_central(skew,shape,lambda,var,c3,c4,m2,m3,m4)
    if(var<=0.0_dp)then;value=ieee_value(0.0_dp,ieee_quiet_nan);else;value=c3/var**1.5_dp;end if
  end function gh_skewness

  real(dp) function gh_excess_kurtosis(skew,shape,lambda) result(value)
    real(dp),intent(in)::skew,shape,lambda
    real(dp)::m2,m3,m4,var,c3,c4
    call gh_raw_central(skew,shape,lambda,var,c3,c4,m2,m3,m4)
    if(var<=0.0_dp)then;value=ieee_value(0.0_dp,ieee_quiet_nan);else;value=c4/(var*var)-3.0_dp;end if
  end function gh_excess_kurtosis

  subroutine gh_raw_central(skew,shape,lambda,var,c3,c4,raw2,raw3,raw4)
    real(dp),intent(in)::skew,shape,lambda;real(dp),intent(out)::var,c3,c4,raw2,raw3,raw4
    real(dp)::a,b,d,m0,chi,psi,w1,w2,w3,w4,r1
    integer::status
    call paramgh(skew,shape,lambda,a,b,d,m0,status)
    if(status/=0)then;var=-1.0_dp;c3=0.0_dp;c4=0.0_dp;raw2=0.0_dp;raw3=0.0_dp;raw4=0.0_dp;return;end if
    chi=d*d;psi=max(a*a-b*b,0.0_dp)
    w1=gig_raw_moment(1.0_dp,lambda,chi,psi);w2=gig_raw_moment(2.0_dp,lambda,chi,psi)
    w3=gig_raw_moment(3.0_dp,lambda,chi,psi);w4=gig_raw_moment(4.0_dp,lambda,chi,psi)
    r1=b*w1;raw2=w1+b*b*w2;raw3=3.0_dp*b*w2+b**3*w3;raw4=3.0_dp*w2+6.0_dp*b*b*w3+b**4*w4
    var=raw2-r1*r1;c3=raw3-3.0_dp*r1*raw2+2.0_dp*r1**3;c4=raw4-4.0_dp*r1*raw3+6.0_dp*r1*r1*raw2-3.0_dp*r1**4
  end subroutine gh_raw_central

  real(dp) function ghst_skewness(skew,shape) result(value)
    real(dp),intent(in)::skew,shape
    real(dp)::b,d,m,b2,d2
    integer::status
    if(shape<=6.0_dp)then;value=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    call paramghst(skew,shape,b,d,m,status);b2=b*b;d2=d*d
    value=2.0_dp*sqrt(shape-4.0_dp)*b*d/(2.0_dp*b2*d2+(shape-2.0_dp)*(shape-4.0_dp))**1.5_dp* &
      (3.0_dp*(shape-2.0_dp)+8.0_dp*b2*d2/(shape-6.0_dp))
  end function ghst_skewness

  real(dp) function ghst_excess_kurtosis(skew,shape) result(value)
    real(dp),intent(in)::skew,shape
    real(dp)::b,d,m,b2,d2,k1
    integer::status
    if(shape<=8.0_dp)then;value=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    call paramghst(skew,shape,b,d,m,status);b2=b*b;d2=d*d;k1=6.0_dp/(2.0_dp*b2*d2+(shape-2.0_dp)*(shape-4.0_dp))**2
    value=k1*((shape-2.0_dp)**2*(shape-4.0_dp)+16.0_dp*b2*d2*(shape-2.0_dp)*(shape-4.0_dp)/(shape-6.0_dp)+ &
      8.0_dp*b2*b2*d2*d2*(5.0_dp*shape-22.0_dp)/((shape-6.0_dp)*(shape-8.0_dp)))
  end function ghst_excess_kurtosis

  function authorized_domain(distribution,max_kurt,n,lambda) result(result)
    character(len=*),intent(in)::distribution
    real(dp),intent(in),optional::max_kurt,lambda
    integer,intent(in),optional::n
    type(authorized_domain_result)::result
    real(dp)::kmax,lam,target,best_error,sval,hval,skv,kv
    real(dp)::slow,sup,hlow,hup
    integer::nn,i,j,k,ns,nh
    kmax=30.0_dp;if(present(max_kurt))kmax=max_kurt
    lam=1.0_dp;if(present(lambda))lam=lambda
    nn=25;if(present(n))nn=max(2,n)
    select case(distribution_id(distribution))
    case(dist_sstd);slow=1.0_dp;sup=12.0_dp;hlow=4.01_dp;hup=100.0_dp
    case(dist_nig);slow=0.0_dp;sup=0.98_dp;hlow=0.05_dp;hup=30.0_dp
    case(dist_gh);slow=0.0_dp;sup=0.98_dp;hlow=0.25_dp;hup=30.0_dp
    case(dist_jsu);slow=0.0_dp;sup=15.0_dp;hlow=0.15_dp;hup=20.0_dp
    case(dist_ghst);slow=0.0_dp;sup=50.0_dp;hlow=8.01_dp;hup=60.0_dp
    case default;result%status=tsd_invalid_argument;result%message='authorized domain is defined for sstd, nig, gh, jsu, or ghst';return
    end select
    allocate(result%skewness(nn),result%kurtosis(nn),result%skew_parameter(nn),result%shape_parameter(nn))
    ns=64;nh=96
    do i=1,nn
      target=3.0_dp+(kmax-3.0_dp)*real(i-1,dp)/real(nn-1,dp)
      best_error=huge(1.0_dp);result%skewness(i)=0.0_dp;result%kurtosis(i)=target
      do j=0,ns
        sval=slow+(sup-slow)*real(j,dp)/real(ns,dp)
        do k=0,nh
          hval=hlow+(hup-hlow)*real(k,dp)/real(nh,dp)
          skv=dskewness(distribution,sval,hval,lam);kv=3.0_dp+dkurtosis(distribution,sval,hval,lam)
          if(.not.ieee_is_finite(skv).or..not.ieee_is_finite(kv))cycle
          if(abs(kv-target)<best_error .or. (abs(kv-target)<=best_error*1.02_dp .and. skv>result%skewness(i)))then
            best_error=abs(kv-target);result%skewness(i)=skv;result%kurtosis(i)=kv;result%skew_parameter(i)=sval;result%shape_parameter(i)=hval
          end if
        end do
      end do
    end do
    result%status=tsd_success;result%message='grid-constrained numerical boundary'
  end function authorized_domain

end module tsd_moments
