module mev_distributions
  use mev_kinds, only: dp, pi
  use mev_math, only: normal_cdf, normal_quantile, reg_gamma_p, gamma_quantile, &
      reg_beta, beta_quantile, expm1_mev
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
  implicit none
  private
  public :: dgev, pgev, qgev, rgev
  public :: dgp, pgp, qgp, rgp
  public :: degp, pegp, qegp, regp
  public :: egp_g_cdf, egp_g_pdf, egp_g_quantile
  public :: gev_endpoint_lower, gev_endpoint_upper, gp_endpoint_upper
contains

  pure real(dp) function gev_endpoint_lower(loc,scale,shape) result(v)
    real(dp),intent(in)::loc,scale,shape
    if(shape>0.0_dp) then
      v=loc-scale/shape
    else
      v=-ieee_value(0.0_dp,ieee_positive_inf)
    end if
  end function

  pure real(dp) function gev_endpoint_upper(loc,scale,shape) result(v)
    real(dp),intent(in)::loc,scale,shape
    if(shape<0.0_dp) then
      v=loc-scale/shape
    else
      v=ieee_value(0.0_dp,ieee_positive_inf)
    end if
  end function

  pure real(dp) function gp_endpoint_upper(loc,scale,shape) result(v)
    real(dp),intent(in)::loc,scale,shape
    if(shape<0.0_dp) then
      v=loc-scale/shape
    else
      v=ieee_value(0.0_dp,ieee_positive_inf)
    end if
  end function

  pure real(dp) function dgev(x,loc,scale,shape,log_density) result(v)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::loc,scale,shape
    logical,intent(in),optional::log_density
    real(dp)::mu,sig,xi,z,t,lv
    logical::lg
    mu=0.0_dp;if(present(loc))mu=loc
    sig=1.0_dp;if(present(scale))sig=scale
    xi=0.0_dp;if(present(shape))xi=shape
    lg=.false.;if(present(log_density))lg=log_density
    if(sig<=0.0_dp)then
      v=ieee_value(0.0_dp,ieee_quiet_nan);return
    end if
    z=(x-mu)/sig
    if(abs(xi)<1.0e-10_dp)then
      lv=-log(sig)-z-exp(-z)
    else
      t=1.0_dp+xi*z
      if(t<=0.0_dp)then
        lv=-ieee_value(0.0_dp,ieee_positive_inf)
      else
        lv=-log(sig)-t**(-1.0_dp/xi)-(1.0_dp/xi+1.0_dp)*log(t)
      end if
    end if
    if(lg)then;v=lv;else;v=exp(lv);end if
  end function dgev

  pure real(dp) function pgev(q,loc,scale,shape,lower_tail,log_p) result(v)
    real(dp),intent(in)::q
    real(dp),intent(in),optional::loc,scale,shape
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::mu,sig,xi,z,t,lp,p
    logical::lt,lg
    mu=0.0_dp;if(present(loc))mu=loc
    sig=1.0_dp;if(present(scale))sig=scale
    xi=0.0_dp;if(present(shape))xi=shape
    lt=.true.;if(present(lower_tail))lt=lower_tail
    lg=.false.;if(present(log_p))lg=log_p
    if(sig<=0.0_dp)then
      v=ieee_value(0.0_dp,ieee_quiet_nan);return
    end if
    z=(q-mu)/sig
    if(abs(xi)>1.0e-6_dp)then
      t=1.0_dp+xi*z
      if(t<=0.0_dp)then
        if(xi>0.0_dp)then;lp=-ieee_value(0.0_dp,ieee_positive_inf);else;lp=0.0_dp;end if
      else
        lp=-t**(-1.0_dp/xi)
      end if
    else
      lp=-exp(-z+0.5_dp*xi*z*z)
    end if
    if(lt)then
      if(lg)then;v=lp;else;v=exp(lp);end if
    else
      p=-expm1_mev(lp)
      if(lg)then;v=log(p);else;v=p;end if
    end if
  end function pgev

  pure real(dp) function qgev(prob,loc,scale,shape,lower_tail,log_p) result(v)
    real(dp),intent(in)::prob
    real(dp),intent(in),optional::loc,scale,shape
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::mu,sig,xi,p,lp,lml
    logical::lt,lg
    mu=0.0_dp;if(present(loc))mu=loc
    sig=1.0_dp;if(present(scale))sig=scale
    xi=0.0_dp;if(present(shape))xi=shape
    lt=.true.;if(present(lower_tail))lt=lower_tail
    lg=.false.;if(present(log_p))lg=log_p
    if(sig<=0.0_dp)then;v=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    if(lg)then
      if(prob>0.0_dp)then;v=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
      p=exp(prob)
    else
      if(prob<0.0_dp.or.prob>1.0_dp)then;v=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
      p=prob
    end if
    if(.not.lt)p=1.0_dp-p
    if(p<=0.0_dp)then;v=gev_endpoint_lower(mu,sig,xi);return;end if
    if(p>=1.0_dp)then;v=gev_endpoint_upper(mu,sig,xi);return;end if
    lp=log(p);lml=log(-lp)
    if(abs(xi)<=1.0e-10_dp)then
      v=mu-sig*lml
    else if(abs(xi)<=1.0e-6_dp)then
      v=mu-sig*lml*(1.0_dp-0.5_dp*lml*xi)
    else
      v=mu+sig*((-lp)**(-xi)-1.0_dp)/xi
    end if
  end function qgev

  subroutine rgev(n,x,loc,scale,shape)
    integer,intent(in)::n
    real(dp),intent(out)::x(n)
    real(dp),intent(in),optional::loc,scale,shape
    real(dp)::u
    integer::i
    do i=1,n
      call random_number(u)
      x(i)=qgev(u,loc,scale,shape)
    end do
  end subroutine rgev

  pure real(dp) function dgp(x,loc,scale,shape,log_density) result(v)
    real(dp),intent(in)::x
    real(dp),intent(in),optional::loc,scale,shape
    logical,intent(in),optional::log_density
    real(dp)::mu,sig,xi,z,t,lv
    logical::lg
    mu=0.0_dp;if(present(loc))mu=loc
    sig=1.0_dp;if(present(scale))sig=scale
    xi=0.0_dp;if(present(shape))xi=shape
    lg=.false.;if(present(log_density))lg=log_density
    if(sig<=0.0_dp)then;v=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    z=(x-mu)/sig
    if(z<0.0_dp)then
      lv=-ieee_value(0.0_dp,ieee_positive_inf)
    else if(abs(xi)<1.0e-10_dp)then
      lv=-log(sig)-z
    else
      t=1.0_dp+xi*z
      if(t<=0.0_dp)then
        lv=-ieee_value(0.0_dp,ieee_positive_inf)
      else
        lv=-log(sig)-(1.0_dp/xi+1.0_dp)*log(t)
      end if
    end if
    if(lg)then;v=lv;else;v=exp(lv);end if
  end function dgp

  pure real(dp) function pgp(q,loc,scale,shape,lower_tail,log_p) result(v)
    real(dp),intent(in)::q
    real(dp),intent(in),optional::loc,scale,shape
    logical,intent(in),optional::lower_tail,log_p
    real(dp)::mu,sig,xi,z,t,p
    logical::lt,lg
    mu=0.0_dp;if(present(loc))mu=loc
    sig=1.0_dp;if(present(scale))sig=scale
    xi=0.0_dp;if(present(shape))xi=shape
    lt=.true.;if(present(lower_tail))lt=lower_tail
    lg=.false.;if(present(log_p))lg=log_p
    if(sig<=0.0_dp)then;v=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    z=max(q-mu,0.0_dp)/sig
    if(abs(xi)<1.0e-8_dp)then
      p=1.0_dp-exp(-z)
    else
      t=1.0_dp+xi*z
      if(t<=0.0_dp.and.xi<0.0_dp)then
        p=1.0_dp
      else
        p=1.0_dp-exp((-1.0_dp/xi)*log(max(t,tiny(1.0_dp))))
      end if
    end if
    if(q<mu)p=0.0_dp
    if(.not.lt)p=1.0_dp-p
    if(lg)then;v=log(p);else;v=p;end if
  end function pgp

  pure real(dp) function qgp(prob,loc,scale,shape,lower_tail) result(v)
    real(dp),intent(in)::prob
    real(dp),intent(in),optional::loc,scale,shape
    logical,intent(in),optional::lower_tail
    real(dp)::mu,sig,xi,p,surv
    logical::lt
    mu=0.0_dp;if(present(loc))mu=loc
    sig=1.0_dp;if(present(scale))sig=scale
    xi=0.0_dp;if(present(shape))xi=shape
    lt=.true.;if(present(lower_tail))lt=lower_tail
    if(prob<0.0_dp.or.prob>1.0_dp.or.sig<=0.0_dp)then
      v=ieee_value(0.0_dp,ieee_quiet_nan);return
    end if
    p=prob
    if(lt)then;surv=1.0_dp-p;else;surv=p;end if
    if(surv>=1.0_dp)then;v=mu;return;end if
    if(surv<=0.0_dp)then;v=gp_endpoint_upper(mu,sig,xi);return;end if
    if(abs(xi)<1.0e-10_dp)then
      v=mu-sig*log(surv)
    else
      v=mu+sig*(surv**(-xi)-1.0_dp)/xi
    end if
  end function qgp

  subroutine rgp(n,x,loc,scale,shape)
    integer,intent(in)::n
    real(dp),intent(out)::x(n)
    real(dp),intent(in),optional::loc,scale,shape
    real(dp)::u
    integer::i
    do i=1,n
      call random_number(u)
      x(i)=qgp(u,loc,scale,shape)
    end do
  end subroutine rgp

  pure real(dp) function beta_logpdf(x,a,b) result(lv)
    real(dp),intent(in)::x,a,b
    if(x<=0.0_dp.or.x>=1.0_dp.or.a<=0.0_dp.or.b<=0.0_dp)then
      lv=-ieee_value(0.0_dp,ieee_positive_inf)
    else
      lv=(a-1.0_dp)*log(x)+(b-1.0_dp)*log(1.0_dp-x) &
         -log_gamma(a)-log_gamma(b)+log_gamma(a+b)
    end if
  end function beta_logpdf

  pure real(dp) function gamma_logpdf(x,a) result(lv)
    real(dp),intent(in)::x,a
    if(x<=0.0_dp.or.a<=0.0_dp)then
      lv=-ieee_value(0.0_dp,ieee_positive_inf)
    else
      lv=(a-1.0_dp)*log(x)-x-log_gamma(a)
    end if
  end function gamma_logpdf

  pure real(dp) function logistic_cdf(x,loc,scale) result(p)
    real(dp),intent(in)::x,loc,scale
    real(dp)::z
    z=(x-loc)/scale
    if(z>=0.0_dp)then;p=1.0_dp/(1.0_dp+exp(-z));else;p=exp(z)/(1.0_dp+exp(z));end if
  end function

  pure real(dp) function logistic_quantile(p,loc,scale) result(x)
    real(dp),intent(in)::p,loc,scale
    x=loc+scale*log(p/(1.0_dp-p))
  end function

  pure real(dp) function egp_g_cdf(x,kappa,shape,model) result(v)
    real(dp),intent(in)::x,kappa,shape
    character(len=*),intent(in)::model
    real(dp)::u,a,cst
    u=max(0.0_dp,min(1.0_dp,x))
    select case(trim(model))
    case('pt-beta')
      if(kappa<=0.0_dp.or.abs(shape)<1.0e-15_dp)then
        v=ieee_value(0.0_dp,ieee_quiet_nan)
      else
        v=reg_beta(1.0_dp-(1.0_dp-u)**abs(shape),kappa,1.0_dp/abs(shape))
      end if
    case('pt-gamma')
      if(kappa<=0.0_dp)then;v=ieee_value(0.0_dp,ieee_quiet_nan);else
        if(u>=1.0_dp)then;v=1.0_dp;else;v=reg_gamma_p(kappa,-log(1.0_dp-u));end if
      end if
    case('pt-power')
      v=u**kappa
    case('gj-beta')
      a=1.0_dp/32.0_dp
      cst=reg_beta(0.5_dp,kappa,kappa)-reg_beta(a,kappa,kappa)
      v=(reg_beta((0.5_dp-a)*u+a,kappa,kappa)-reg_beta(a,kappa,kappa))/cst
      v=max(0.0_dp,min(1.0_dp,v))
    case('gj-tnorm')
      if(kappa<=0.0_dp)then
        v=u
      else
        cst=0.5_dp-normal_cdf(-sqrt(kappa))
        v=(normal_cdf((u-1.0_dp)*sqrt(kappa))-normal_cdf(-sqrt(kappa)))/cst
      end if
    case('exptilt')
      if(abs(kappa-1.0_dp)<1.0e-12_dp)then;v=u;else;v=(kappa**u-1.0_dp)/(kappa-1.0_dp);end if
    case('logist')
      if(kappa<=1.0e-15_dp)then
        v=u
      else
        cst=logistic_cdf(0.0_dp,1.0_dp,1.0_dp/kappa)
        v=(logistic_cdf(u,1.0_dp,1.0_dp/kappa)-cst)/(0.5_dp-cst)
        v=max(0.0_dp,min(1.0_dp,v))
      end if
    case default
      v=ieee_value(0.0_dp,ieee_quiet_nan)
    end select
  end function egp_g_cdf

  pure real(dp) function egp_g_pdf(x,kappa,shape,model,log_density) result(v)
    real(dp),intent(in)::x,kappa,shape
    character(len=*),intent(in)::model
    logical,intent(in),optional::log_density
    real(dp)::lv,a,cst,z
    logical::lg
    lg=.false.;if(present(log_density))lg=log_density
    if(x<0.0_dp.or.x>1.0_dp)then
      lv=-ieee_value(0.0_dp,ieee_positive_inf)
    else
      select case(trim(model))
      case('pt-beta')
        if(abs(kappa-1.0_dp)<1.0e-12_dp)then
          lv=0.0_dp
        else if(x<=0.0_dp.or.x>=1.0_dp.or.abs(shape)<1.0e-15_dp)then
          lv=-ieee_value(0.0_dp,ieee_positive_inf)
        else
          z=1.0_dp-(1.0_dp-x)**abs(shape)
          lv=log(abs(shape))+(abs(shape)-1.0_dp)*log(1.0_dp-x) &
             +beta_logpdf(z,kappa,1.0_dp/abs(shape))
        end if
      case('pt-gamma')
        if(x>=1.0_dp)then;lv=-ieee_value(0.0_dp,ieee_positive_inf);else
          z=-log(max(1.0_dp-x,tiny(1.0_dp)))
          if(z<=0.0_dp.and.kappa>1.0_dp)then
            lv=-ieee_value(0.0_dp,ieee_positive_inf)
          else
            lv=gamma_logpdf(max(z,tiny(1.0_dp)),kappa)-log(max(1.0_dp-x,tiny(1.0_dp)))
          end if
        end if
      case('pt-power')
        if(abs(kappa-1.0_dp)<1.0e-12_dp)then;lv=0.0_dp;else
          if(x<=0.0_dp)then;lv=-ieee_value(0.0_dp,ieee_positive_inf);else;lv=log(kappa)+(kappa-1.0_dp)*log(x);end if
        end if
      case('gj-beta')
        a=1.0_dp/32.0_dp
        cst=reg_beta(0.5_dp,kappa,kappa)-reg_beta(a,kappa,kappa)
        lv=log(0.5_dp-a)+beta_logpdf((0.5_dp-a)*x+a,kappa,kappa)-log(cst)
      case('gj-tnorm')
        if(kappa<=1.0e-15_dp)then;lv=0.0_dp;else
          cst=0.5_dp-normal_cdf(-sqrt(kappa))
          lv=0.5_dp*log(kappa)-0.5_dp*log(2.0_dp*pi) &
             -0.5_dp*kappa*(x-1.0_dp)**2-log(cst)
        end if
      case('exptilt')
        if(abs(kappa-1.0_dp)<1.0e-12_dp)then;lv=0.0_dp;else
          lv=x*log(kappa)+log(abs(log(kappa)))-log(abs(kappa-1.0_dp))
        end if
      case('logist')
        if(kappa<=1.0e-15_dp)then;lv=0.0_dp;else
          cst=0.5_dp-logistic_cdf(0.0_dp,1.0_dp,1.0_dp/kappa)
          z=(x-1.0_dp)*kappa
          lv=log(kappa)-z-2.0_dp*log(1.0_dp+exp(-z))-log(cst)
        end if
      case default
        lv=ieee_value(0.0_dp,ieee_quiet_nan)
      end select
    end if
    if(lg)then;v=lv;else;v=exp(lv);end if
  end function egp_g_pdf

  pure real(dp) function egp_g_quantile(p,kappa,shape,model) result(v)
    real(dp),intent(in)::p,kappa,shape
    character(len=*),intent(in)::model
    real(dp)::u,a,cst,z
    u=max(0.0_dp,min(1.0_dp,p))
    select case(trim(model))
    case('pt-beta')
      if(abs(shape)<1.0e-15_dp)then;v=ieee_value(0.0_dp,ieee_quiet_nan);else
        z=beta_quantile(u,kappa,1.0_dp/abs(shape))
        v=1.0_dp-(1.0_dp-z)**(1.0_dp/abs(shape))
      end if
    case('pt-gamma')
      v=1.0_dp-exp(-gamma_quantile(u,kappa))
    case('pt-power')
      v=u**(1.0_dp/kappa)
    case('gj-beta')
      a=1.0_dp/32.0_dp
      cst=reg_beta(0.5_dp,kappa,kappa)-reg_beta(a,kappa,kappa)
      z=beta_quantile(u*cst+reg_beta(a,kappa,kappa),kappa,kappa)
      v=(z-a)/(0.5_dp-a)
    case('gj-tnorm')
      if(kappa<=1.0e-15_dp)then;v=u;else
        cst=0.5_dp-normal_cdf(-sqrt(kappa))
        z=u*cst+normal_cdf(-sqrt(kappa))
        v=1.0_dp+normal_quantile(z)/sqrt(kappa)
      end if
    case('exptilt')
      if(abs(kappa-1.0_dp)<1.0e-12_dp)then;v=u;else;v=log(1.0_dp+u*(kappa-1.0_dp))/log(kappa);end if
    case('logist')
      if(kappa<=1.0e-15_dp)then;v=u;else
        cst=logistic_cdf(0.0_dp,1.0_dp,1.0_dp/kappa)
        v=logistic_quantile(u*(0.5_dp-cst)+cst,1.0_dp,1.0_dp/kappa)
      end if
    case default
      v=ieee_value(0.0_dp,ieee_quiet_nan)
    end select
    v=max(0.0_dp,min(1.0_dp,v))
  end function egp_g_quantile

  pure real(dp) function pegp(q,scale,shape,kappa,model,lower_tail,log_p) result(v)
    real(dp),intent(in)::q,scale,shape,kappa
    character(len=*),intent(in)::model
    logical,intent(in),optional::lower_tail,log_p
    logical::lt,lg
    real(dp)::p
    lt=.true.;if(present(lower_tail))lt=lower_tail
    lg=.false.;if(present(log_p))lg=log_p
    if(scale<=0.0_dp.or.kappa<=0.0_dp)then;v=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    if((trim(model)=='pt-beta'.or.trim(model)=='pt-gamma').and.abs(shape)<1.0e-8_dp)then
      if(q<=0.0_dp)then;p=0.0_dp;else;p=reg_gamma_p(kappa,q/scale);end if
    else
      p=egp_g_cdf(pgp(q,0.0_dp,scale,shape),kappa,shape,model)
    end if
    if(.not.lt)p=1.0_dp-p
    if(lg)then;v=log(p);else;v=p;end if
  end function pegp

  pure real(dp) function degp(x,scale,shape,kappa,model,log_density) result(v)
    real(dp),intent(in)::x,scale,shape,kappa
    character(len=*),intent(in)::model
    logical,intent(in),optional::log_density
    logical::lg
    real(dp)::lv,p
    lg=.false.;if(present(log_density))lg=log_density
    if(scale<=0.0_dp.or.kappa<0.0_dp)then;v=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    if((trim(model)=='pt-beta'.or.trim(model)=='pt-gamma').and.abs(shape)<1.0e-8_dp)then
      if(x<=0.0_dp)then;lv=-ieee_value(0.0_dp,ieee_positive_inf);else
        lv=(kappa-1.0_dp)*log(x/scale)-x/scale-log_gamma(kappa)-log(scale)
      end if
    else
      p=pgp(x,0.0_dp,scale,shape)
      lv=egp_g_pdf(p,kappa,shape,model,.true.)+dgp(x,0.0_dp,scale,shape,.true.)
    end if
    if(lg)then;v=lv;else;v=exp(lv);end if
  end function degp

  pure real(dp) function qegp(prob,scale,shape,kappa,model,lower_tail,log_p) result(v)
    real(dp),intent(in)::prob,scale,shape,kappa
    character(len=*),intent(in)::model
    logical,intent(in),optional::lower_tail,log_p
    logical::lt,lg
    real(dp)::p,u
    lt=.true.;if(present(lower_tail))lt=lower_tail
    lg=.false.;if(present(log_p))lg=log_p
    p=prob;if(lg)p=exp(p)
    if(.not.lt)p=1.0_dp-p
    if(p<0.0_dp.or.p>1.0_dp.or.scale<=0.0_dp.or.kappa<0.0_dp)then
      v=ieee_value(0.0_dp,ieee_quiet_nan);return
    end if
    if((trim(model)=='pt-beta'.or.trim(model)=='pt-gamma').and.abs(shape)<1.0e-8_dp)then
      v=gamma_quantile(p,kappa,scale)
    else
      u=egp_g_quantile(p,kappa,shape,model)
      v=qgp(u,0.0_dp,scale,shape)
    end if
  end function qegp

  subroutine regp(n,x,scale,shape,kappa,model)
    integer,intent(in)::n
    real(dp),intent(out)::x(n)
    real(dp),intent(in)::scale,shape,kappa
    character(len=*),intent(in)::model
    real(dp)::u
    integer::i
    do i=1,n
      call random_number(u)
      x(i)=qegp(u,scale,shape,kappa,model)
    end do
  end subroutine regp

end module mev_distributions
