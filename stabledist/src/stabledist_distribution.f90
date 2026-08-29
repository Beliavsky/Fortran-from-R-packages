! Computational translation of R package stabledist 0.7-2.
! Original R implementation: Diethelm Wuertz; accuracy/fixes: Martin Maechler.
! Stable density/probability integrals follow J. P. Nolan's formulas.
! SPDX-License-Identifier: GPL-2.0-or-later
module stabledist_distribution
   use r_compat, only: dp, normal_cdf, dnorm, qnorm, dcauchy, pcauchy, qcauchy, &
      r_gamma, r_lgamma, runif_vec
   use stabledist_numerics, only: bisect_root, golden_max, integrate_split
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, &
      ieee_positive_inf, ieee_negative_inf
   implicit none
   private

   real(dp), parameter :: pi = acos(-1.0_dp)
   real(dp), parameter :: pi2 = 0.5_dp*pi
   real(dp), parameter :: log_huge_dp = log(huge(1.0_dp))
   real(dp), parameter :: log_tiny_dp = log(tiny(1.0_dp))
   real(dp), parameter :: large_exp_arg = -log(tiny(1.0_dp))
   real(dp), parameter :: alpha_small_dstable = 1.0e-17_dp

   public :: dstable, pstable, qstable, rstable, rstable_varying, stable_mode
   public :: c_stable_tail, dpareto, ppareto, stable_omega

   interface dstable
      module procedure dstable_scalar, dstable_vec
   end interface
   interface pstable
      module procedure pstable_scalar, pstable_vec
   end interface
   interface qstable
      module procedure qstable_scalar, qstable_vec
   end interface
   interface c_stable_tail
      module procedure c_stable_tail_scalar, c_stable_tail_vec
   end interface
   interface dpareto
      module procedure dpareto_scalar, dpareto_vec
   end interface
   interface ppareto
      module procedure ppareto_scalar, ppareto_vec
   end interface

contains

   pure function nan_dp() result(x)
      real(dp) :: x
      x = ieee_value(0.0_dp, ieee_quiet_nan)
   end function nan_dp

   pure function pos_inf() result(x)
      real(dp) :: x
      x = ieee_value(0.0_dp, ieee_positive_inf)
   end function pos_inf

   pure function neg_inf() result(x)
      real(dp) :: x
      x = ieee_value(0.0_dp, ieee_negative_inf)
   end function neg_inf

   pure function log1p_s(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      if (abs(x) < 1.0e-8_dp) then
         y = x*(1.0_dp + x*(-0.5_dp + x*(1.0_dp/3.0_dp + x*(-0.25_dp + 0.2_dp*x))))
      else
         y = log(1.0_dp+x)
      end if
   end function log1p_s

   pure function expm1_s(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      if (abs(x) < 1.0e-6_dp) then
         y = x*(1.0_dp + x*(0.5_dp + x*(1.0_dp/6.0_dp + x*(1.0_dp/24.0_dp + x/120.0_dp))))
      else
         y = exp(x)-1.0_dp
      end if
   end function expm1_s

   pure function tanpi2_scalar(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      real(dp) :: xr
      xr = anint(x)
      if (x /= 0.0_dp .and. x == xr) then
         if (modulo(nint(x),4) == 1) then
            y = pos_inf()
         else if (modulo(nint(x),4) == 3) then
            y = neg_inf()
         else
            y = tan(pi2*x)
         end if
      else
         y = tan(pi2*x)
      end if
   end function tanpi2_scalar

   pure function cospi2_scalar(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      if (x == anint(x)) then
         if (x == 0.0_dp) then
            y = 1.0_dp
         else
            y = 0.0_dp
         end if
      else
         y = cos(pi2*x)
      end if
   end function cospi2_scalar

   pure function stable_omega(gamma, alpha) result(om)
      real(dp), intent(in) :: gamma, alpha
      real(dp) :: om
      if (alpha /= anint(alpha)) then
         om = tan(pi2*alpha)
      else if (alpha == 1.0_dp) then
         if (gamma > 0.0_dp) then
            om = (2.0_dp/pi)*log(gamma)
         else
            om = neg_inf()
         end if
      else
         om = 0.0_dp
      end if
   end function stable_omega

   pure function c_stable_tail_scalar(alpha, log_) result(ans)
      real(dp), intent(in) :: alpha
      logical, intent(in), optional :: log_
      real(dp) :: ans
      logical :: lg
      lg = .false.
      if (present(log_)) lg = log_
      if (alpha < 0.0_dp .or. alpha > 2.0_dp) then
         ans = nan_dp()
         return
      end if
      if (alpha == 0.0_dp) then
         ans = merge(-log(2.0_dp), 0.5_dp, lg)
      else if (alpha == 2.0_dp) then
         ans = merge(neg_inf(), 0.0_dp, lg)
      else if (lg) then
         ans = r_lgamma(alpha)-log(pi)+log(sin(alpha*pi2))
      else
         ans = r_gamma(alpha)/pi*sin(alpha*pi2)
      end if
   end function c_stable_tail_scalar

   pure function c_stable_tail_vec(alpha, log_) result(ans)
      real(dp), intent(in) :: alpha(:)
      logical, intent(in), optional :: log_
      real(dp), allocatable :: ans(:)
      integer :: i
      allocate(ans(size(alpha)))
      do i=1,size(alpha)
         ans(i)=c_stable_tail_scalar(alpha(i),log_)
      end do
   end function c_stable_tail_vec

   pure function dpareto_scalar(x, alpha, beta, log_) result(ans)
      real(dp), intent(in) :: x, alpha, beta
      logical, intent(in), optional :: log_
      real(dp) :: ans, xx, bb
      logical :: lg
      lg=.false.
      if(present(log_)) lg=log_
      xx=x
      bb=beta
      if (xx < 0.0_dp) then
         xx=-xx
         bb=-bb
      end if
      if (xx <= 0.0_dp .or. alpha <= 0.0_dp .or. alpha > 2.0_dp .or. abs(bb)>1.0_dp) then
         ans=merge(neg_inf(),0.0_dp,lg)
         return
      end if
      if (1.0_dp+bb <= 0.0_dp .or. c_stable_tail_scalar(alpha) == 0.0_dp) then
         ans=merge(neg_inf(),0.0_dp,lg)
         return
      end if
      if (lg) then
         ans=log(alpha)+log1p_s(bb)+c_stable_tail_scalar(alpha,.true.)-(1.0_dp+alpha)*log(xx)
      else
         ans=alpha*(1.0_dp+bb)*c_stable_tail_scalar(alpha)*xx**(-(1.0_dp+alpha))
      end if
   end function dpareto_scalar

   pure function dpareto_vec(x, alpha, beta, log_) result(ans)
      real(dp), intent(in) :: x(:), alpha, beta
      logical, intent(in), optional :: log_
      real(dp), allocatable :: ans(:)
      integer :: i
      allocate(ans(size(x)))
      do i=1,size(x)
         ans(i)=dpareto_scalar(x(i),alpha,beta,log_)
      end do
   end function dpareto_vec

   pure function ppareto_scalar(x, alpha, beta, lower_tail, log_p) result(ans)
      real(dp), intent(in) :: x, alpha, beta
      logical, intent(in), optional :: lower_tail, log_p
      real(dp) :: ans, ifail
      logical :: lower, lp
      lower=.true.
      if(present(lower_tail)) lower=lower_tail
      lp=.false.
      if(present(log_p)) lp=log_p
      if (x < 0.0_dp) then
         ans=nan_dp()
         return ! upstream marks negative-x implementation FIXME/incorrect
      end if
      if (x == 0.0_dp) then
         ans=merge(merge(neg_inf(),0.0_dp,lp),merge(0.0_dp,1.0_dp,lp),lower)
         return
      end if
      ifail=(1.0_dp+beta)*c_stable_tail_scalar(alpha)*x**(-alpha)
      ifail=max(0.0_dp,min(1.0_dp,ifail))
      if (lower) then
         if (lp) then
            ans=log1p_s(-ifail)
         else
            ans=1.0_dp-ifail
         end if
      else
         if (lp) then
            if (ifail==0.0_dp) then
            ans=neg_inf()
            else
            ans=log(ifail)
            end if
         else
            ans=ifail
         end if
      end if
   end function ppareto_scalar

   pure function ppareto_vec(x, alpha, beta, lower_tail, log_p) result(ans)
      real(dp), intent(in) :: x(:), alpha, beta
      logical, intent(in), optional :: lower_tail, log_p
      real(dp), allocatable :: ans(:)
      integer :: i
      allocate(ans(size(x)))
      do i=1,size(x)
         ans(i)=ppareto_scalar(x(i),alpha,beta,lower_tail,log_p)
      end do
   end function ppareto_vec

   pure function levy_s0_density(z,beta,lg) result(f)
      real(dp),intent(in)::z,beta
      logical,intent(in)::lg
      real(dp)::f,y
      if(beta>0.0_dp)then
      y=z+1.0_dp
      else
      y=1.0_dp-z
      end if
      if(y<=0.0_dp)then
         f=merge(neg_inf(),0.0_dp,lg)
      else if(lg)then
         f=-0.5_dp*log(2.0_dp*pi)-0.5_dp/y-1.5_dp*log(y)
      else
         f=exp(-0.5_dp/y)/(sqrt(2.0_dp*pi)*y**1.5_dp)
      end if
   end function levy_s0_density

   pure function levy_s0_cdf(z,beta) result(f)
      real(dp),intent(in)::z,beta
      real(dp)::f,y,fl
      if(beta>0.0_dp)then
      y=z+1.0_dp
      else
      y=1.0_dp-z
      end if
      if(y<=0.0_dp)then
         fl=0.0_dp
      else
         fl=erfc(sqrt(0.5_dp/y))
      end if
      if(beta>0.0_dp)then
      f=fl
      else
      f=1.0_dp-fl
      end if
   end function levy_s0_cdf

   recursive subroutine pm_transform(alpha,beta,gamma,delta,pm,gam,del)
      real(dp), intent(in) :: alpha,beta,gamma,delta
      integer, intent(in) :: pm
      real(dp), intent(out) :: gam,del
      gam=gamma
      del=delta
      if (pm==1) then
         del=delta+beta*gamma*stable_omega(gamma,alpha)
      else if (pm==2) then
         gam=alpha**(-1.0_dp/alpha)*gamma
         del=delta-gam*stable_mode(alpha,beta)
      end if
   end subroutine pm_transform

   recursive function dstable_scalar(x, alpha, beta, gamma, delta, pm, log_, tol, zeta_tol, subdivisions) result(ans)
      real(dp), intent(in) :: x, alpha, beta
      real(dp), intent(in), optional :: gamma, delta, tol, zeta_tol
      integer, intent(in), optional :: pm, subdivisions
      logical, intent(in), optional :: log_
      real(dp) :: ans
      real(dp) :: gam0, del0, gam, del, z, t, ztol, betan, zeta, theta0, f
      integer :: prm, nsub
      logical :: lg

      gam0=1.0_dp
      if(present(gamma)) gam0=gamma
      del0=0.0_dp
      if(present(delta)) del0=delta
      prm=0
      if(present(pm)) prm=pm
      lg=.false.
      if(present(log_)) lg=log_
      t=64.0_dp*epsilon(1.0_dp)
      if(present(tol)) t=tol
      nsub=1000
      if(present(subdivisions)) nsub=subdivisions
      if (alpha<=0.0_dp .or. alpha>2.0_dp .or. abs(beta)>1.0_dp .or. gam0<0.0_dp .or. &
          prm<0 .or. prm>2 .or. t<=0.0_dp .or. nsub<=0) then
         ans=nan_dp()
         return
      end if
      if (gam0==0.0_dp) then
         if (x==del0) then
         ans=pos_inf()
         else
         ans=merge(neg_inf(),0.0_dp,lg)
         end if
         return
      end if
      call pm_transform(alpha,beta,gam0,del0,prm,gam,del)
      z=(x-del)/gam
      if (alpha==2.0_dp) then
         f=dnorm(z,mean=0.0_dp,sd=sqrt(2.0_dp),log_=lg)
      else if (alpha==1.0_dp .and. beta==0.0_dp) then
         f=dcauchy(z,location=0.0_dp,scale=1.0_dp,log_=lg)
      else if (alpha==0.5_dp .and. abs(beta)==1.0_dp) then
         f=levy_s0_density(z,beta,lg)
      else if (alpha/=1.0_dp) then
         betan=beta*tan(pi2*alpha)
         zeta=-betan
         theta0=min(max(-pi2,atan(betan)/alpha),pi2)
         if (present(zeta_tol)) then
            ztol=max(0.0_dp,zeta_tol)
         else if (betan==0.0_dp) then
            ztol=0.4e-15_dp
         else if (1.0_dp-abs(beta)<0.01_dp .or. alpha<0.01_dp) then
            ztol=2.0e-15_dp
         else
            ztol=5.0e-5_dp
         end if
         f=density_alpha_ne1(z,zeta,alpha,beta,theta0,lg,t,ztol,nsub)
      else
         if (z>=0.0_dp) then
            f=density_alpha1(abs(z),beta,lg,t,nsub)
         else
            f=density_alpha1(abs(z),-beta,lg,t,nsub)
         end if
      end if
      if ((.not.lg .and. f==0.0_dp) .or. (lg .and. f==neg_inf())) then
         f=dpareto_scalar(z,alpha,beta,lg)
      end if
      if (lg) then
         ans=f-log(gam)
      else
         ans=f/gam
      end if
   end function dstable_scalar

   function dstable_vec(x, alpha, beta, gamma, delta, pm, log_, tol, zeta_tol, subdivisions) result(ans)
      real(dp), intent(in) :: x(:), alpha, beta
      real(dp), intent(in), optional :: gamma(:), delta(:), tol, zeta_tol
      integer, intent(in), optional :: pm, subdivisions
      logical, intent(in), optional :: log_
      real(dp), allocatable :: ans(:)
      real(dp) :: gg,dd
      integer :: i
      allocate(ans(size(x)))
      do i=1,size(x)
         gg=1.0_dp
         dd=0.0_dp
         if(present(gamma)) gg=gamma(1+mod(i-1,size(gamma)))
         if(present(delta)) dd=delta(1+mod(i-1,size(delta)))
         ans(i)=dstable_scalar(x(i),alpha,beta,gg,dd,pm,log_,tol,zeta_tol,subdivisions)
      end do
   end function dstable_vec

   function density_alpha_ne1(xin,zeta,alpha,betain,theta0in,lg,tol,ztol,nsub) result(ans)
      real(dp), intent(in) :: xin,zeta,alpha,betain,theta0in,tol,ztol
      integer, intent(in) :: nsub
      logical, intent(in) :: lg
      real(dp) :: ans
      real(dp) :: x,beta,theta0,xmz,fz,a1,at0,cat0,c2,a,b,theta2
      real(dp) :: ga,gb,gpa,gpb,gpeak,epsp,th1,th3,r1,r2,r3,r4,sumr
      logical :: ok1,ok3,okroot,small_alpha

      x=xin
      beta=betain
      theta0=theta0in
      xmz=abs(x-zeta)
      if (lg) then
         fz=r_lgamma(1.0_dp+1.0_dp/alpha)+log(cos(theta0))- &
            (log(pi)+log1p_s(zeta*zeta)/(2.0_dp*alpha))
      else
         fz=r_gamma(1.0_dp+1.0_dp/alpha)*cos(theta0)/(pi*(1.0_dp+zeta*zeta)**(1.0_dp/(2.0_dp*alpha)))
      end if
      if (ieee_is_finite(x) .and. xmz <= ztol*(ztol+max(abs(x),abs(zeta)))) then
         ans=fz
         return
      end if
      small_alpha=alpha<alpha_small_dstable
      if (x<zeta) then
         theta0=-theta0
         if (small_alpha) then
            beta=-beta
            x=-x
         end if
      end if
      if (small_alpha) then
         if (x<=0.0_dp) then
            ans=merge(neg_inf(),0.0_dp,lg)
         else if (lg) then
            ans=log(alpha)+log1p_s(beta)-log(2.0_dp*x+pi*alpha*beta)
         else
            ans=alpha*(1.0_dp+beta)/(2.0_dp*x+pi*alpha*beta)
         end if
         return
      end if
      a1=alpha-1.0_dp
      at0=alpha*theta0
      cat0=cos(at0)
      c2=alpha/(pi*abs(a1)*xmz)
      a=-theta0
      b=pi2
      ga=stable_g(a+1.0e-12_dp*max(1.0_dp,abs(a)),a1,at0,cat0,xmz,alpha)
      gb=stable_g(b-1.0e-12_dp,a1,at0,cat0,xmz,alpha)
      gpa=g_times_exp_minus_g(ga)
      gpb=g_times_exp_minus_g(gb)
      if ((ga-1.0_dp)*(gb-1.0_dp)<=0.0_dp .and. ieee_is_finite(ga) .and. ieee_is_finite(gb)) then
         theta2=bisect_root(root_g,a+1.0e-12_dp*max(1.0_dp,abs(a)),b-1.0e-12_dp,epsilon(1.0_dp),200,okroot)
         if (.not.okroot) theta2=0.5_dp*(a+b)
      else
         if (gpa>=gpb) then
         theta2=a
         else
         theta2=b
         end if
      end if
      gpeak=g1(theta2)
      epsp=1.0e-4_dp
      ok1=.false.
      ok3=.false.
      th1=theta2
      th3=theta2
      if (theta2>a .and. gpeak>epsp .and. g1(a)<epsp) then
         th1=bisect_root(root_left,a,theta2,max(tol,epsilon(1.0_dp)),200,ok1)
      end if
      if (theta2<b .and. gpeak>epsp .and. g1(b)<epsp) then
         th3=bisect_root(root_right,theta2,b,max(tol,epsilon(1.0_dp)),200,ok3)
      end if
      r1=0.0_dp
      r2=0.0_dp
      r3=0.0_dp
      r4=0.0_dp
      if (ok1) then
         r1=integrate_split(g1,a,th1,tol,nsub,8)
         r2=integrate_split(g1,th1,theta2,tol,nsub,8)
      else if (theta2>a) then
         r2=integrate_split(g1,a,theta2,tol,nsub,12)
      end if
      if (ok3) then
         r3=integrate_split(g1,theta2,th3,tol,nsub,8)
         r4=integrate_split(g1,th3,b,tol,nsub,8)
      else if (b>theta2) then
         r3=integrate_split(g1,theta2,b,tol,nsub,12)
      end if
      sumr=max(0.0_dp,r1+r2+r3+r4)
      if (sumr==0.0_dp) then
         ans=merge(neg_inf(),0.0_dp,lg)
      else if (lg) then
         ans=log(c2)+log(sumr)
      else
         ans=c2*sumr
      end if
   contains
      function root_g(th) result(v)
         real(dp),intent(in)::th
         real(dp)::v
         v=stable_g(th,a1,at0,cat0,xmz,alpha)-1.0_dp
      end function root_g
      function g1(th) result(v)
         real(dp),intent(in)::th
         real(dp)::v
         v=g_times_exp_minus_g(stable_g(th,a1,at0,cat0,xmz,alpha))
      end function g1
      function root_left(th) result(v)
         real(dp),intent(in)::th
         real(dp)::v
         v=g1(th)-epsp
      end function root_left
      function root_right(th) result(v)
         real(dp),intent(in)::th
         real(dp)::v
         v=g1(th)-epsp
      end function root_right
   end function density_alpha_ne1

   pure function stable_g(th,a1,at0,cat0,xmz,alpha) result(g)
      real(dp),intent(in)::th,a1,at0,cat0,xmz,alpha
      real(dp)::g,att,ct,sa,ca,logg
      if (abs(pi2-sign(1.0_dp,a1)*th)<64.0_dp*epsilon(1.0_dp)) then
         g=0.0_dp
         return
      end if
      att=at0+alpha*th
      ct=cos(th)
      sa=sin(att)
      ca=cos(att-th)
      if (ct<=0.0_dp .or. cat0<=0.0_dp .or. ca<=0.0_dp) then
         g=0.0_dp
         return
      end if
      if (sa<=0.0_dp) then
         if (a1>0.0_dp) then
         g=pos_inf()
         else
         g=0.0_dp
         end if
         return
      end if
      logg=(log(cat0)+log(ct)+alpha*(log(xmz)-log(sa)))/a1+log(ca)
      if (logg>log_huge_dp) then
         g=pos_inf()
      else if (logg<log_tiny_dp) then
         g=0.0_dp
      else
         g=exp(logg)
      end if
   end function stable_g

   pure function g_times_exp_minus_g(g) result(v)
      real(dp),intent(in)::g
      real(dp)::v
      if (.not.ieee_is_finite(g) .or. g>large_exp_arg .or. g<=0.0_dp) then
         v=0.0_dp
      else
         v=g*exp(-g)
      end if
   end function g_times_exp_minus_g

   function density_alpha1(x,beta,lg,tol,nsub) result(ans)
      real(dp),intent(in)::x,beta,tol
      integer,intent(in)::nsub
      logical,intent(in)::lg
      real(dp)::ans,i2b,p2b,ea,u2,r1,r2
      logical::ok
      if (beta==0.0_dp) then
         ans=dcauchy(x,location=0.0_dp,scale=1.0_dp,log_=lg)
         return
      end if
      i2b=1.0_dp/(2.0_dp*beta)
      p2b=pi*i2b
      ea=-p2b*x
      if (.not.ieee_is_finite(ea)) then
         ans=merge(neg_inf(),0.0_dp,lg)
         return
      end if
      u2=bisect_root(rootg,-1.0_dp+1.0e-12_dp,1.0_dp-1.0e-12_dp,max(tol,epsilon(1.0_dp)),200,ok)
      if(.not.ok) u2=0.0_dp
      r1=integrate_split(g2,-1.0_dp,u2,tol,nsub,16)
      r2=integrate_split(g2,u2,1.0_dp,tol,nsub,16)
      if(r1+r2<=0.0_dp) then
         ans=merge(neg_inf(),0.0_dp,lg)
      else if(lg) then
         ans=log(pi2)+log(abs(i2b))+log(r1+r2)
      else
         ans=pi2*abs(i2b)*(r1+r2)
      end if
   contains
      function rootg(u) result(v)
         real(dp),intent(in)::u
         real(dp)::v
         v=g_alpha1(u,beta,ea)-1.0_dp
      end function rootg
      function g2(u) result(v)
         real(dp),intent(in)::u
         real(dp)::v
         v=g_times_exp_minus_g(g_alpha1(u,beta,ea))
      end function g2
   end function density_alpha1

   pure function g_alpha1(u,beta,ea) result(g)
      real(dp),intent(in)::u,beta,ea
      real(dp)::g,i2b,p2b,th,h,ratio,c,tt,logg,u0
      u0=-sign(1.0_dp,beta)
      if(abs(u-u0)<1.0e-10_dp) then
      g=0.0_dp
      return
      end if
      if(abs(u)>=1.0_dp) then
         if (u==u0) then
         g=0.0_dp
         else
         g=pos_inf()
         end if
         return
      end if
      i2b=1.0_dp/(2.0_dp*beta)
      p2b=pi*i2b
      th=u*pi2
      h=p2b+th
      ratio=h/p2b
      c=cospi2_scalar(u)
      tt=tanpi2_scalar(u)
      if(ratio<=0.0_dp .or. c<=0.0_dp) then
      g=0.0_dp
      return
      end if
      logg=log(ratio)+ea+h*tt-log(c)
      if(logg>log_huge_dp) then
      g=pos_inf()
      else if(logg<log_tiny_dp) then
      g=0.0_dp
      else
      g=exp(logg)
      end if
   end function g_alpha1

   recursive function pstable_scalar(q,alpha,beta,gamma,delta,pm,lower_tail,log_p,tol,subdivisions) result(ans)
      real(dp),intent(in)::q,alpha,beta
      real(dp),intent(in),optional::gamma,delta,tol
      integer,intent(in),optional::pm,subdivisions
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::ans,gam0,del0,gam,del,x,t,zeta,theta0,f1,r,bb
      integer::prm,nsub
      logical::lower,lp,fin_supp,use_lower,give_i,use_l
      gam0=1.0_dp
      if(present(gamma))gam0=gamma
      del0=0.0_dp
      if(present(delta))del0=delta
      prm=0
      if(present(pm))prm=pm
      lower=.true.
      if(present(lower_tail))lower=lower_tail
      lp=.false.
      if(present(log_p))lp=log_p
      t=64.0_dp*epsilon(1.0_dp)
      if(present(tol))t=tol
      nsub=1000
      if(present(subdivisions))nsub=subdivisions
      if(alpha<=0.0_dp.or.alpha>2.0_dp.or.abs(beta)>1.0_dp.or.gam0<0.0_dp.or.prm<0.or.prm>2.or.t<=0.or.nsub<=0)then
         ans=nan_dp()
         return
      end if
      if(gam0==0.0_dp)then
         r=merge(1.0_dp,0.0_dp,q>=del0)
         ans=prob_output(r,lower,lp)
         return
      end if
      call pm_transform(alpha,beta,gam0,del0,prm,gam,del)
      x=(q-del)/gam
      if(alpha==2.0_dp)then
         r=normal_cdf(x/sqrt(2.0_dp))
         ans=prob_output(r,lower,lp)
         return
      else if(alpha==1.0_dp.and.beta==0.0_dp)then
         r=pcauchy(x,location=0.0_dp,scale=1.0_dp)
         ans=prob_output(r,lower,lp)
         return
      else if(alpha==0.5_dp.and.abs(beta)==1.0_dp)then
         r=levy_s0_cdf(x,beta)
         ans=prob_output(r,lower,lp)
         return
      end if
      if(alpha/=1.0_dp)then
         zeta=-beta*tan(pi2*alpha)
         theta0=min(max(-pi2,atan(-zeta)/alpha),pi2)
         fin_supp=(abs(beta)==1.0_dp.and.alpha<1.0_dp)
         if(fin_supp)then
            if(beta==1.0_dp.and.x<=zeta)then
            ans=prob_output(0.0_dp,lower,lp)
            return
            end if
            if(beta==-1.0_dp.and.x>=zeta)then
            ans=prob_output(1.0_dp,lower,lp)
            return
            end if
         end if
         if(abs(x-zeta)<2.0_dp*epsilon(1.0_dp))then
            r=0.5_dp-theta0/pi
            ans=prob_output(r,lower,lp)
            return
         end if
         use_lower=((x>zeta.and.lower).or.(x<zeta.and..not.lower))
         give_i=(.not.use_lower.and.alpha>1.0_dp)
         f1=prob_alpha_ne1(x,zeta,alpha,theta0,give_i,t,nsub)
         if(give_i)then
            if(lp)then
            if(f1<=0.0_dp)then
            ans=neg_inf()
            else
            ans=log(f1)
            end if
            else
            ans=f1
            end if
         else
            ans=prob_output(f1,use_lower,lp)
         end if
      else
         bb=beta
         if(bb>=0.0_dp)then
            use_l=lower
         else
            bb=-bb
            x=-x
            use_l=.not.lower
         end if
         give_i=(.not.use_l.and..not.lp)
         if(give_i)use_l=.true.
         f1=prob_alpha1(x,bb,t,nsub,give_i)
         ans=prob_output(f1,use_l,lp)
      end if
   end function pstable_scalar

   function pstable_vec(q,alpha,beta,gamma,delta,pm,lower_tail,log_p,tol,subdivisions) result(ans)
      real(dp),intent(in)::q(:),alpha,beta
      real(dp),intent(in),optional::gamma(:),delta(:),tol
      integer,intent(in),optional::pm,subdivisions
      logical,intent(in),optional::lower_tail,log_p
      real(dp),allocatable::ans(:)
      real(dp)::gg,dd
      integer::i
      allocate(ans(size(q)))
      do i=1,size(q)
         gg=1.0_dp
         dd=0.0_dp
         if(present(gamma))gg=gamma(1+mod(i-1,size(gamma)))
         if(present(delta))dd=delta(1+mod(i-1,size(delta)))
         ans(i)=pstable_scalar(q(i),alpha,beta,gg,dd,pm,lower_tail,log_p,tol,subdivisions)
      end do
   end function pstable_vec

   pure function prob_output(f,lower,lp) result(ans)
      real(dp),intent(in)::f
      logical,intent(in)::lower,lp
      real(dp)::ans,p
      p=max(0.0_dp,min(1.0_dp,f))
      if(lower)then
         if(lp)then
         if(p<=0.0_dp)then
         ans=neg_inf()
         else
         ans=log(p)
         end if
         else
         ans=p
         end if
      else
         if(lp)then
            if(p>=1.0_dp)then
            ans=neg_inf()
            else
            ans=log1p_s(-p)
            end if
         else
         ans=1.0_dp-p
         end if
      end if
   end function prob_output

   function prob_alpha_ne1(x,zeta,alpha,theta0in,give_i,tol,nsub) result(ans)
      real(dp),intent(in)::x,zeta,alpha,theta0in,tol
      logical,intent(in)::give_i
      integer,intent(in)::nsub
      real(dp)::ans,theta0,xmz,a1,at0,cat0,a,b,r,c1,c3
      if(.not.ieee_is_finite(x))then
      ans=merge(0.0_dp,1.0_dp,give_i)
      return
      end if
      theta0=theta0in
      xmz=abs(x-zeta)
      if(x<zeta)theta0=-theta0
      a1=alpha-1.0_dp
      at0=alpha*theta0
      cat0=cos(at0)
      a=-theta0+1.0e-12_dp*max(1.0_dp,abs(theta0))
      b=pi2-1.0e-12_dp
      if(a>=b)then
      a=-theta0
      b=pi2
      end if
      r=integrate_split(gexp,a,b,tol,nsub,32)
      if(give_i)then
         ans=max(0.0_dp,min(1.0_dp,r/pi))
      else
         if(alpha<1.0_dp)then
         c1=0.5_dp-theta0/pi
         else
         c1=1.0_dp
         end if
         c3=sign(1.0_dp,1.0_dp-alpha)/pi
         ans=max(0.0_dp,min(1.0_dp,c1+c3*r))
      end if
   contains
      function gexp(th) result(v)
         real(dp),intent(in)::th
         real(dp)::v,g
         g=stable_g(th,a1,at0,cat0,xmz,alpha)
         if(.not.ieee_is_finite(g).or.g>large_exp_arg)then
         v=0.0_dp
         else
         v=exp(-max(0.0_dp,g))
         end if
      end function gexp
   end function prob_alpha_ne1

   function prob_alpha1(x,beta,tol,nsub,give_i) result(ans)
      real(dp),intent(in)::x,beta,tol
      integer,intent(in)::nsub
      logical,intent(in)::give_i
      real(dp)::ans,i2b,p2b,ea,r
      if(beta==0.0_dp)then
         ans=pcauchy(x,location=0.0_dp,scale=1.0_dp)
         return
      end if
      i2b=1.0_dp/(2.0_dp*beta)
      p2b=pi*i2b
      ea=-p2b*x
      if(.not.ieee_is_finite(ea))then
         if(ea<0.0_dp)then
         ans=merge(0.0_dp,1.0_dp,give_i)
         else
         ans=0.0_dp
         end if
         return
      end if
      r=0.5_dp*integrate_split(gprob,-1.0_dp,1.0_dp,tol,nsub,32)
      if(give_i)then
      ans=-r
      else
      ans=r
      end if
      ans=max(0.0_dp,min(1.0_dp,ans))
   contains
      function gprob(u) result(v)
         real(dp),intent(in)::u
         real(dp)::v,g
         g=g_alpha1(u,beta,ea)
         if(give_i)then
            if(.not.ieee_is_finite(g))then
            v=-1.0_dp
            else
            v=expm1_s(-max(0.0_dp,g))
            end if
         else
            if(.not.ieee_is_finite(g).or.g>large_exp_arg)then
            v=0.0_dp
            else
            v=exp(-max(0.0_dp,g))
            end if
         end if
      end function gprob
   end function prob_alpha1

   recursive function qstable_scalar(p,alpha,beta,gamma,delta,pm,lower_tail,log_p,tol,maxiter,integ_tol,subdivisions) result(ans)
      real(dp),intent(in)::p,alpha,beta
      real(dp),intent(in),optional::gamma,delta,tol,integ_tol
      integer,intent(in),optional::pm,maxiter,subdivisions
      logical,intent(in),optional::lower_tail,log_p
      real(dp)::ans,gam0,del0,gam,del,ptarget,qtol,it,lo,hi,flo,fhi,mid,fm,zeta
      integer::prm,mi,nsub,k
      logical::lower,lp
      gam0=1.0_dp
      if(present(gamma))gam0=gamma
      del0=0.0_dp
      if(present(delta))del0=delta
      prm=0
      if(present(pm))prm=pm
      lower=.true.
      if(present(lower_tail))lower=lower_tail
      lp=.false.
      if(present(log_p))lp=log_p
      qtol=epsilon(1.0_dp)**0.25_dp
      if(present(tol))qtol=tol
      it=1.0e-7_dp
      if(present(integ_tol))it=integ_tol
      mi=1000
      if(present(maxiter))mi=maxiter
      nsub=200
      if(present(subdivisions))nsub=subdivisions
      if(alpha<=0.0_dp.or.alpha>2.0_dp.or.abs(beta)>1.0_dp.or.gam0<0.0_dp.or.prm<0.or.prm>2)then
      ans=nan_dp()
      return
      end if
      if(lp)then
         if(p>0.0_dp)then
         ans=nan_dp()
         return
         end if
         if(lower)then
         ptarget=exp(p)
         else
         ptarget=-expm1_s(p)
         end if
      else
         if(p<0.0_dp.or.p>1.0_dp)then
         ans=nan_dp()
         return
         end if
         ptarget=merge(p,1.0_dp-p,lower)
      end if
      if(ptarget<=0.0_dp)then
      ans=neg_inf()
      return
      end if
      if(ptarget>=1.0_dp)then
      ans=pos_inf()
      return
      end if
      if(gam0==0.0_dp)then
      ans=del0
      return
      end if
      call pm_transform(alpha,beta,gam0,del0,prm,gam,del)
      if(alpha==2.0_dp)then
         ans=qnorm(ptarget,mean=0.0_dp,sd=sqrt(2.0_dp),lower_tail=.true.)*gam+del
         return
      else if(alpha==1.0_dp.and.beta==0.0_dp)then
         ans=qcauchy(ptarget,location=0.0_dp,scale=1.0_dp)*gam+del
         return
      end if
      lo=-1.0_dp
      hi=1.0_dp
      if(alpha<1.0_dp.and.abs(beta)==1.0_dp)then
         zeta=-beta*tan(pi2*alpha)
         if(beta==1.0_dp)lo=zeta
         if(beta==-1.0_dp)hi=zeta
      end if
      flo=pstable_scalar(lo,alpha,beta,1.0_dp,0.0_dp,0,.true.,.false.,it,nsub)-ptarget
      fhi=pstable_scalar(hi,alpha,beta,1.0_dp,0.0_dp,0,.true.,.false.,it,nsub)-ptarget
      do k=1,200
         if(flo<=0.0_dp.and.fhi>=0.0_dp)exit
         if (flo > 0.0_dp) then
            hi=lo
            fhi=flo
            lo=2.0_dp*lo-1.0_dp
            flo=pstable_scalar(lo,alpha,beta,tol=it,subdivisions=nsub)-ptarget
         end if
         if (fhi < 0.0_dp) then
            lo=hi
            flo=fhi
            hi=2.0_dp*hi+1.0_dp
            fhi=pstable_scalar(hi,alpha,beta,tol=it,subdivisions=nsub)-ptarget
         end if
      end do
      if(.not.(flo<=0.0_dp.and.fhi>=0.0_dp))then
      ans=nan_dp()
      return
      end if
      do k=1,mi
         mid=0.5_dp*(lo+hi)
         fm=pstable_scalar(mid,alpha,beta,tol=it,subdivisions=nsub)-ptarget
         if(abs(hi-lo)<=qtol*max(1.0_dp,abs(mid)))exit
         if(fm<0.0_dp)then
         lo=mid
         else
         hi=mid
         end if
      end do
      ans=0.5_dp*(lo+hi)*gam+del
   end function qstable_scalar

   function qstable_vec(p,alpha,beta,gamma,delta,pm,lower_tail,log_p,tol,maxiter,integ_tol,subdivisions) result(ans)
      real(dp),intent(in)::p(:),alpha,beta
      real(dp),intent(in),optional::gamma(:),delta(:),tol,integ_tol
      integer,intent(in),optional::pm,maxiter,subdivisions
      logical,intent(in),optional::lower_tail,log_p
      real(dp),allocatable::ans(:)
      real(dp)::gg,dd
      integer::i
      allocate(ans(size(p)))
      do i=1,size(p)
         gg=1.0_dp
         dd=0.0_dp
         if(present(gamma))gg=gamma(1+mod(i-1,size(gamma)))
         if(present(delta))dd=delta(1+mod(i-1,size(delta)))
         ans(i)=qstable_scalar(p(i),alpha,beta,gg,dd,pm,lower_tail,log_p,tol,maxiter,integ_tol,subdivisions)
      end do
   end function qstable_vec

   function rstable(n,alpha,beta,gamma,delta,pm) result(ans)
      integer,intent(in)::n
      real(dp),intent(in)::alpha,beta
      real(dp),intent(in),optional::gamma,delta
      integer,intent(in),optional::pm
      real(dp),allocatable::ans(:)
      real(dp),allocatable::u1(:),u2(:)
      real(dp)::gam0,del0,gam,del,theta,w,btan,theta0,c,ath,r
      integer::prm,i
      allocate(ans(max(0,n)))
      if(n<=0)return
      gam0=1.0_dp
      if(present(gamma))gam0=gamma
      del0=0.0_dp
      if(present(delta))del0=delta
      prm=0
      if(present(pm))prm=pm
      if(alpha<=0.0_dp.or.alpha>2.0_dp.or.abs(beta)>1.0_dp.or.gam0<0.0_dp.or.prm<0.or.prm>2)then
      ans=nan_dp()
      return
      end if
      call pm_transform(alpha,beta,gam0,del0,prm,gam,del)
      u1=runif_vec(n)
      u2=runif_vec(n)
      do i=1,n
         theta=pi*(u1(i)-0.5_dp)
         w=-log(max(u2(i),tiny(1.0_dp)))
         if(alpha==1.0_dp.and.beta==0.0_dp)then
            r=tan(theta)
         else if(alpha==1.0_dp)then
            ! Stable limit of the Chambers-Mallows-Stuck generator at alpha=1.
            r=(2.0_dp/pi)*((pi2+beta*theta)*tan(theta)- &
               beta*log((pi2*w*cos(theta))/(pi2+beta*theta)))
         else
            btan=beta*tan(pi2*alpha)
            theta0=min(max(-pi2,atan(btan)/alpha),pi2)
            c=(1.0_dp+btan*btan)**(1.0_dp/(2.0_dp*alpha))
            ath=alpha*(theta+theta0)
            r=(c*sin(ath)/(cos(theta)**(1.0_dp/alpha)))* &
               (cos(theta-ath)/w)**((1.0_dp-alpha)/alpha)-btan
         end if
         ans(i)=r*gam+del
      end do
   end function rstable

   function rstable_varying(n,alpha,beta,gamma,delta,pm) result(ans)
      integer,intent(in)::n
      real(dp),intent(in)::alpha,beta,gamma(:),delta(:)
      integer,intent(in),optional::pm
      real(dp),allocatable::ans(:)
      real(dp),allocatable::u1(:),u2(:)
      real(dp)::gg,dd,theta,w,btan,theta0,c,ath,r,fac,sm
      integer::prm,i
      allocate(ans(max(0,n)))
      if(n<=0)return
      prm=0
      if(present(pm))prm=pm
      if(alpha<=0.0_dp.or.alpha>2.0_dp.or.abs(beta)>1.0_dp.or.prm<0.or.prm>2)then
         ans=nan_dp()
         return
      end if
      if(any(gamma<0.0_dp).or.size(gamma)==0.or.size(delta)==0)then
         ans=nan_dp()
         return
      end if
      fac=alpha**(-1.0_dp/alpha)
      sm=0.0_dp
      if(prm==2)sm=stable_mode(alpha,beta)
      u1=runif_vec(n)
      u2=runif_vec(n)
      do i=1,n
         gg=gamma(1+mod(i-1,size(gamma)))
         dd=delta(1+mod(i-1,size(delta)))
         if(prm==1)then
            dd=dd+beta*gg*stable_omega(gg,alpha)
         else if(prm==2)then
            gg=fac*gg
            dd=dd-gg*sm
         end if
         theta=pi*(u1(i)-0.5_dp)
         w=-log(max(u2(i),tiny(1.0_dp)))
         if(alpha==1.0_dp.and.beta==0.0_dp)then
            r=tan(theta)
         else if(alpha==1.0_dp)then
            r=(2.0_dp/pi)*((pi2+beta*theta)*tan(theta)- &
               beta*log((pi2*w*cos(theta))/(pi2+beta*theta)))
         else
            btan=beta*tan(pi2*alpha)
            theta0=min(max(-pi2,atan(btan)/alpha),pi2)
            c=(1.0_dp+btan*btan)**(1.0_dp/(2.0_dp*alpha))
            ath=alpha*(theta+theta0)
            r=(c*sin(ath)/(cos(theta)**(1.0_dp/alpha)))* &
               (cos(theta-ath)/w)**((1.0_dp-alpha)/alpha)-btan
         end if
         ans(i)=r*gg+dd
      end do
   end function rstable_varying

   function stable_mode(alpha,beta,beta_max,tol) result(mode)
      real(dp),intent(in)::alpha,beta
      real(dp),intent(in),optional::beta_max,tol
      real(dp)::mode,bm,t,bb,a,b
      if(alpha<=0.0_dp.or.alpha>2.0_dp.or.abs(beta)>1.0_dp)then
      mode=nan_dp()
      return
      end if
      if(alpha*beta==0.0_dp.or.alpha==2.0_dp)then
      mode=0.0_dp
      return
      end if
      bm=1.0_dp-1.0e-11_dp
      if(present(beta_max))bm=beta_max
      t=epsilon(1.0_dp)**0.25_dp
      if(present(tol))t=tol
      bb=beta
      if(bb>bm)bb=bm
      a=-0.7_dp*sign(1.0_dp,bb)
      b=0.0_dp
      mode=golden_max(dens,a,b,t,200)
   contains
      function dens(x) result(v)
         real(dp),intent(in)::x
         real(dp)::v
         v=dstable_scalar(x,alpha,bb,pm=0,tol=max(64.0_dp*epsilon(1.0_dp),t),subdivisions=1000)
      end function dens
   end function stable_mode

end module stabledist_distribution
