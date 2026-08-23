module skewhyperbolic_distribution
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
  use skewhyp_math, only: dp, pi, bessel_k, gamma_rand, normal_rand, adaptive_simpson, &
       binom_coeff, double_factorial_odd
  implicit none
  private
  public :: dskewhyp, dskewhyp_vec, ddskewhyp, pskewhyp, pskewhyp_vec, qskewhyp, qskewhyp_vec
  public :: rskewhyp, skewhyp_mean, skewhyp_var, skewhyp_skew, skewhyp_kurt, skewhyp_mode
  public :: skewhyp_moment, skewhyp_calc_range, skewhyp_check_pars
contains
  pure logical function skewhyp_check_pars(mu,delta,beta,nu)
    real(dp), intent(in) :: mu,delta,beta,nu
    skewhyp_check_pars = (delta > 0.0_dp .and. nu > 0.0_dp .and. &
                          .not. ieee_is_nan(mu) .and. .not. ieee_is_nan(beta))
  end function skewhyp_check_pars

  function log_dskewhyp(x,mu,delta,beta,nu) result(ld)
    real(dp), intent(in) :: x,mu,delta,beta,nu
    real(dp) :: ld, r, z, order, kv, ab
    if (.not. skewhyp_check_pars(mu,delta,beta,nu)) then
      ld = ieee_value(0.0_dp,ieee_quiet_nan)
      return
    end if
    r = sqrt(delta*delta + (x-mu)**2)
    if (abs(beta) <= sqrt(epsilon(1.0_dp))) then
      ld = log_gamma(0.5_dp*(nu+1.0_dp)) - 0.5_dp*log(pi) - log(delta) - &
           log_gamma(0.5_dp*nu) - 0.5_dp*(nu+1.0_dp)*log(1.0_dp + ((x-mu)/delta)**2)
    else
      ab = abs(beta)
      z = ab*r
      order = 0.5_dp*(nu+1.0_dp)
      kv = bessel_k(order,z)
      if (kv <= 0.0_dp) then
        ld = -huge(1.0_dp)
        return
      end if
      ld = 0.5_dp*(1.0_dp-nu)*log(2.0_dp) + nu*log(delta) + &
           order*log(ab) + log(kv) + beta*(x-mu) - log_gamma(0.5_dp*nu) - &
           0.5_dp*log(pi) - order*log(r)
    end if
  end function log_dskewhyp

  function dskewhyp(x,mu,delta,beta,nu,log_density) result(d)
    real(dp), intent(in) :: x,mu,delta,beta,nu
    logical, intent(in), optional :: log_density
    real(dp) :: d, ld
    logical :: lg
    lg = .false.; if (present(log_density)) lg = log_density
    ld = log_dskewhyp(x,mu,delta,beta,nu)
    if (lg) then
      d = ld
    else if (ld < log(tiny(1.0_dp))) then
      d = 0.0_dp
    else
      d = exp(ld)
    end if
  end function dskewhyp

  subroutine dskewhyp_vec(x,mu,delta,beta,nu,d,log_density)
    real(dp), intent(in) :: x(:),mu,delta,beta,nu
    real(dp), intent(out) :: d(size(x))
    logical, intent(in), optional :: log_density
    integer :: i
    do i=1,size(x)
      d(i)=dskewhyp(x(i),mu,delta,beta,nu,log_density)
    end do
  end subroutine dskewhyp_vec

  function ddskewhyp(x,mu,delta,beta,nu) result(dd)
    real(dp), intent(in) :: x,mu,delta,beta,nu
    real(dp) :: dd,h
    h = epsilon(1.0_dp)**(1.0_dp/3.0_dp)*max(1.0_dp,abs(x),delta)
    dd = (dskewhyp(x+h,mu,delta,beta,nu)-dskewhyp(x-h,mu,delta,beta,nu))/(2.0_dp*h)
  end function ddskewhyp

  function pskewhyp(q,mu,delta,beta,nu,lower_tail,log_p,tol) result(p)
    real(dp), intent(in) :: q,mu,delta,beta,nu
    logical, intent(in), optional :: lower_tail,log_p
    real(dp), intent(in), optional :: tol
    real(dp) :: p, mode, reltol, tail, q0
    logical :: lower, lp
    lower=.true.; lp=.false.; reltol=1.0e-8_dp
    if(present(lower_tail)) lower=lower_tail
    if(present(log_p)) lp=log_p
    if(present(tol)) reltol=tol
    if(q <= -huge(1.0_dp)/2) then
      p=0.0_dp
    else if(q >= huge(1.0_dp)/2) then
      p=1.0_dp
    else
      mode=skewhyp_mode(mu,delta,beta,nu)
      q0=q
      if(q <= mode) then
        p = adaptive_simpson(left_fun,0.0_dp,1.0_dp,reltol)
      else
        tail = adaptive_simpson(right_fun,0.0_dp,1.0_dp,reltol)
        p = 1.0_dp-tail
      end if
      p=max(0.0_dp,min(1.0_dp,p))
    end if
    if(.not.lower) p=1.0_dp-p
    if(lp) then
      if(p<=0.0_dp) then; p=-huge(1.0_dp); else; p=log(p); end if
    end if
  contains
    function left_fun(u) result(y)
      real(dp),intent(in)::u
      real(dp)::y,x,jac
      if(u<=0.0_dp) then; y=0.0_dp; return; endif
      x=q0-(1.0_dp-u)/u
      jac=1.0_dp/(u*u)
      y=dskewhyp(x,mu,delta,beta,nu)*jac
    end function left_fun
    function right_fun(u) result(y)
      real(dp),intent(in)::u
      real(dp)::y,x,jac
      if(u<=0.0_dp) then; y=0.0_dp; return; endif
      x=q0+(1.0_dp-u)/u
      jac=1.0_dp/(u*u)
      y=dskewhyp(x,mu,delta,beta,nu)*jac
    end function right_fun
  end function pskewhyp

  subroutine pskewhyp_vec(q,mu,delta,beta,nu,p,lower_tail,log_p,tol)
    real(dp),intent(in)::q(:),mu,delta,beta,nu
    real(dp),intent(out)::p(size(q))
    logical,intent(in),optional::lower_tail,log_p
    real(dp),intent(in),optional::tol
    integer::i
    do i=1,size(q)
      p(i)=pskewhyp(q(i),mu,delta,beta,nu,lower_tail,log_p,tol)
    enddo
  end subroutine pskewhyp_vec

  function qskewhyp(prob,mu,delta,beta,nu,lower_tail,log_p,tol) result(q)
    real(dp),intent(in)::prob,mu,delta,beta,nu
    logical,intent(in),optional::lower_tail,log_p
    real(dp),intent(in),optional::tol
    real(dp)::q,p,lo,hi,mid,ptol,mode,step
    logical::lower,lp
    integer::i
    lower=.true.;lp=.false.;ptol=1.0e-8_dp
    if(present(lower_tail))lower=lower_tail
    if(present(log_p))lp=log_p
    if(present(tol))ptol=tol
    p=prob
    if(lp)p=exp(p)
    if(.not.lower)p=1.0_dp-p
    if(p<=0.0_dp)then;q=-huge(1.0_dp);return;endif
    if(p>=1.0_dp)then;q=huge(1.0_dp);return;endif
    mode=skewhyp_mode(mu,delta,beta,nu)
    step=max(delta,1.0_dp)
    lo=mode-step;hi=mode+step
    do while(pskewhyp(lo,mu,delta,beta,nu,tol=ptol)>p)
      step=2.0_dp*step;lo=mode-step
      if(step>1.0e8_dp)exit
    enddo
    step=max(delta,1.0_dp)
    do while(pskewhyp(hi,mu,delta,beta,nu,tol=ptol)<p)
      step=2.0_dp*step;hi=mode+step
      if(step>1.0e8_dp)exit
    enddo
    do i=1,70
      mid=0.5_dp*(lo+hi)
      if(pskewhyp(mid,mu,delta,beta,nu,tol=ptol)<p)then;lo=mid;else;hi=mid;endif
    enddo
    q=0.5_dp*(lo+hi)
  end function qskewhyp

  subroutine qskewhyp_vec(p,mu,delta,beta,nu,q,lower_tail,log_p,tol)
    real(dp),intent(in)::p(:),mu,delta,beta,nu
    real(dp),intent(out)::q(size(p))
    logical,intent(in),optional::lower_tail,log_p
    real(dp),intent(in),optional::tol
    integer::i
    do i=1,size(p)
      q(i)=qskewhyp(p(i),mu,delta,beta,nu,lower_tail,log_p,tol)
    enddo
  end subroutine qskewhyp_vec

  subroutine rskewhyp(x,mu,delta,beta,nu)
    real(dp),intent(out)::x(:)
    real(dp),intent(in)::mu,delta,beta,nu
    real(dp)::g,w
    integer::i
    if(.not.skewhyp_check_pars(mu,delta,beta,nu))then
      x=ieee_value(0.0_dp,ieee_quiet_nan);return
    endif
    do i=1,size(x)
      g=gamma_rand(0.5_dp*nu)*(2.0_dp/(delta*delta))
      w=1.0_dp/g
      x(i)=mu+beta*w+sqrt(w)*normal_rand()
    enddo
  end subroutine rskewhyp

  function invgamma_moment(s,delta,nu) result(m)
    integer,intent(in)::s
    real(dp),intent(in)::delta,nu
    real(dp)::m
    if(nu<=2.0_dp*real(s,dp))then
      m=ieee_value(0.0_dp,ieee_quiet_nan)
    else
      m=(0.5_dp*delta*delta)**real(s,dp)*exp(log_gamma(0.5_dp*nu-real(s,dp))-log_gamma(0.5_dp*nu))
    endif
  end function invgamma_moment

  function mu_moment(order,delta,beta,nu) result(m)
    integer,intent(in)::order
    real(dp),intent(in)::delta,beta,nu
    real(dp)::m,term,wm
    integer::r,j,s
    m=0.0_dp
    do r=0,order/2
      j=2*r
      s=order-r
      wm=invgamma_moment(s,delta,nu)
      if(ieee_is_nan(wm))then;m=wm;return;endif
      term=binom_coeff(order,j)*beta**(order-j)*double_factorial_odd(r)*wm
      m=m+term
    enddo
  end function mu_moment

  function skewhyp_moment(order,mu,delta,beta,nu,about) result(m)
    integer,intent(in)::order
    real(dp),intent(in)::mu,delta,beta,nu
    real(dp),intent(in),optional::about
    real(dp)::m,a,mm
    integer::k
    a=0.0_dp;if(present(about))a=about
    m=0.0_dp
    do k=0,order
      if(k==0)then;mm=1.0_dp;else;mm=mu_moment(k,delta,beta,nu);endif
      if(ieee_is_nan(mm))then;m=mm;return;endif
      m=m+binom_coeff(order,k)*(mu-a)**(order-k)*mm
    enddo
  end function skewhyp_moment

  function skewhyp_mean(mu,delta,beta,nu) result(m)
    real(dp),intent(in)::mu,delta,beta,nu;real(dp)::m
    if(nu<=2.0_dp)then;m=ieee_value(0.0_dp,ieee_quiet_nan);else;m=mu+beta*delta**2/(nu-2.0_dp);endif
  end function skewhyp_mean
  function skewhyp_var(mu,delta,beta,nu) result(v)
    real(dp),intent(in)::mu,delta,beta,nu;real(dp)::v
    if(ieee_is_nan(mu))then;v=ieee_value(0.0_dp,ieee_quiet_nan);return;endif
    if(nu<=4.0_dp)then;v=ieee_value(0.0_dp,ieee_quiet_nan);else
      v=2.0_dp*beta**2*delta**4/((nu-2.0_dp)**2*(nu-4.0_dp))+delta**2/(nu-2.0_dp)
    endif
  end function skewhyp_var
  function skewhyp_skew(mu,delta,beta,nu) result(s)
    real(dp),intent(in)::mu,delta,beta,nu;real(dp)::s
    if(ieee_is_nan(mu))then;s=ieee_value(0.0_dp,ieee_quiet_nan);return;endif
    if(nu<=6.0_dp)then;s=ieee_value(0.0_dp,ieee_quiet_nan);else
      s=(2.0_dp*sqrt(nu-4.0_dp)*beta*delta/(2.0_dp*beta**2*delta**2+(nu-2.0_dp)*(nu-4.0_dp))**1.5_dp)* &
        (3.0_dp*(nu-2.0_dp)+8.0_dp*beta**2*delta**2/(nu-6.0_dp))
    endif
  end function skewhyp_skew
  function skewhyp_kurt(mu,delta,beta,nu) result(k)
    real(dp),intent(in)::mu,delta,beta,nu;real(dp)::k
    if(ieee_is_nan(mu))then;k=ieee_value(0.0_dp,ieee_quiet_nan);return;endif
    if(nu<=8.0_dp)then;k=ieee_value(0.0_dp,ieee_quiet_nan);else
      k=6.0_dp/(2.0_dp*beta**2*delta**2+(nu-2.0_dp)*(nu-4.0_dp))**2 * &
        ((nu-2.0_dp)**2*(nu-4.0_dp)+16.0_dp*beta**2*delta**2*(nu-2.0_dp)*(nu-4.0_dp)/(nu-6.0_dp)+ &
        8.0_dp*beta**4*delta**4*(5.0_dp*nu-22.0_dp)/((nu-6.0_dp)*(nu-8.0_dp)))
    endif
  end function skewhyp_kurt

  function skewhyp_mode(mu,delta,beta,nu) result(mode)
    real(dp),intent(in)::mu,delta,beta,nu
    real(dp)::mode,a,b,c,d1,d2
    integer::i
    if(abs(beta)<=sqrt(epsilon(1.0_dp)))then;mode=mu;return;endif
    a=mu-10.0_dp*max(delta,abs(beta)*delta*delta/max(nu,1.0_dp),1.0_dp)
    b=mu+10.0_dp*max(delta,abs(beta)*delta*delta/max(nu,1.0_dp),1.0_dp)
    do i=1,100
      c=a+(b-a)/3.0_dp; mode=b-(b-a)/3.0_dp
      d1=dskewhyp(c,mu,delta,beta,nu,log_density=.true.)
      d2=dskewhyp(mode,mu,delta,beta,nu,log_density=.true.)
      if(d1<d2)then;a=c;else;b=mode;endif
    enddo
    mode=0.5_dp*(a+b)
  end function skewhyp_mode

  subroutine skewhyp_calc_range(mu,delta,beta,nu,range,tol,density)
    real(dp),intent(in)::mu,delta,beta,nu
    real(dp),intent(out)::range(2)
    real(dp),intent(in),optional::tol
    logical,intent(in),optional::density
    real(dp)::tt,mode,lo,hi,mid
    logical::dens
    integer::i
    tt=1.0e-5_dp;if(present(tol))tt=tol
    dens=.true.;if(present(density))dens=density
    mode=skewhyp_mode(mu,delta,beta,nu)
    lo=mode-max(delta,1.0_dp);hi=mode+max(delta,1.0_dp)
    if(dens)then
      do while(dskewhyp(lo,mu,delta,beta,nu)>tt);lo=mode-2.0_dp*(mode-lo);enddo
      do while(dskewhyp(hi,mu,delta,beta,nu)>tt);hi=mode+2.0_dp*(hi-mode);enddo
      call bisect_density(lo,mode,range(1));call bisect_density(mode,hi,range(2))
    else
      do while(pskewhyp(lo,mu,delta,beta,nu)>tt);lo=mode-2.0_dp*(mode-lo);enddo
      do while(1.0_dp-pskewhyp(hi,mu,delta,beta,nu)>tt);hi=mode+2.0_dp*(hi-mode);enddo
      do i=1,60;mid=0.5_dp*(lo+mode);if(pskewhyp(mid,mu,delta,beta,nu)<tt)then;lo=mid;else;mode=mid;endif;enddo
      range(1)=0.5_dp*(lo+mode)
      mode=skewhyp_mode(mu,delta,beta,nu)
      do i=1,60;mid=0.5_dp*(mode+hi);if(1.0_dp-pskewhyp(mid,mu,delta,beta,nu)>tt)then;mode=mid;else;hi=mid;endif;enddo
      range(2)=0.5_dp*(mode+hi)
    endif
  contains
    subroutine bisect_density(a0,b0,r)
      real(dp),intent(in)::a0,b0;real(dp),intent(out)::r
      real(dp)::aa,bb,mm
      integer::j
      aa=a0;bb=b0
      do j=1,60
        mm=0.5_dp*(aa+bb)
        if((dskewhyp(aa,mu,delta,beta,nu)-tt)*(dskewhyp(mm,mu,delta,beta,nu)-tt)<=0.0_dp)then
          bb=mm
        else
          aa=mm
        endif
      enddo
      r=0.5_dp*(aa+bb)
    end subroutine bisect_density
  end subroutine skewhyp_calc_range
end module skewhyperbolic_distribution
