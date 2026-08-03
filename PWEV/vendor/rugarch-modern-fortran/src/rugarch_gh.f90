! Generalized-hyperbolic distribution support for rugarch-modern-fortran.
! Derived from the standardized parameterizations and computational behavior of
! the GPL-3 R package rugarch 1.5-6. This independent Fortran translation is
! distributed under GPL-3.0-only; see LICENSE, NOTICE, and ORIGIN.md.
module rugarch_gh
   use rugarch_kinds, only : dp
   use rugarch_math, only : pi, student_t_pdf
   use rugarch_rng, only : random_normal, random_gamma
   implicit none
   private

   public :: bessel_k_log, bessel_k_scaled_log
   public :: gh_parameters, ghst_parameters
   public :: dgh_raw, dsgh, psgh, qsgh, rsgh
   public :: dsnig, psnig, qsnig, rsnig
   public :: dsghst, psghst, qsghst, rsghst
   public :: random_gig

contains

   pure elemental function log_cosh_stable(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value, ax
      ax = abs(x)
      value = ax - log(2.0_dp) + log(1.0_dp+exp(-2.0_dp*ax))
   end function log_cosh_stable

   function bessel_k_scaled_log(x, nu) result(value)
      real(dp), intent(in) :: x, nu
      real(dp) :: value
      real(dp) :: anu, tmax, t, lf, peak, h, sumv, corr
      integer :: i, n, quiet

      if (x <= 0.0_dp) then
         value = huge(1.0_dp)
         return
      end if
      anu = abs(nu)

      if (x > 50.0_dp) then
         corr = 1.0_dp + (4.0_dp*anu*anu-1.0_dp)/(8.0_dp*x) + &
            (4.0_dp*anu*anu-1.0_dp)*(4.0_dp*anu*anu-9.0_dp)/(128.0_dp*x*x)
         value = 0.5_dp*log(pi/(2.0_dp*x)) + log(max(corr,tiny(1.0_dp)))
         return
      end if

      if (x < 1.0e-5_dp) then
         if (anu > 1.0e-6_dp) then
            value = log(0.5_dp) + log_gamma(anu) + anu*log(2.0_dp/x) + x
         else
            value = log(max(-log(0.5_dp*x)-0.5772156649015329_dp,tiny(1.0_dp))) + x
         end if
         return
      end if

      peak = 0.0_dp
      tmax = 0.0_dp
      quiet = 0
      do i = 1, 400
         tmax = 0.1_dp*real(i,dp)
         lf = -x*(cosh(tmax)-1.0_dp) + log_cosh_stable(anu*tmax)
         peak = max(peak,lf)
         if (lf < peak-48.0_dp .and. tmax > 2.0_dp) then
            quiet = quiet+1
            if (quiet >= 4) exit
         else
            quiet = 0
         end if
      end do
      tmax = max(tmax,4.0_dp)
      n = max(120,2*ceiling(18.0_dp*tmax))
      if (mod(n,2) /= 0) n=n+1
      h = tmax/real(n,dp)

      peak = -huge(1.0_dp)
      do i = 0, n
         t = h*real(i,dp)
         lf = -x*(cosh(t)-1.0_dp) + log_cosh_stable(anu*t)
         peak = max(peak,lf)
      end do

      sumv = 0.0_dp
      do i = 0, n
         t = h*real(i,dp)
         lf = -x*(cosh(t)-1.0_dp) + log_cosh_stable(anu*t)
         if (i==0 .or. i==n) then
            sumv = sumv + exp(lf-peak)
         else if (mod(i,2)==0) then
            sumv = sumv + 2.0_dp*exp(lf-peak)
         else
            sumv = sumv + 4.0_dp*exp(lf-peak)
         end if
      end do
      value = peak + log(h*sumv/3.0_dp)
   end function bessel_k_scaled_log

   function bessel_k_log(x, nu) result(value)
      real(dp), intent(in) :: x, nu
      real(dp) :: value
      value = bessel_k_scaled_log(x,nu)-x
   end function bessel_k_log

   function kappa_gh(x, lambda) result(value)
      real(dp), intent(in) :: x, lambda
      real(dp) :: value
      if (abs(lambda+0.5_dp) < 1.0e-13_dp) then
         value = 1.0_dp/x
      else
         value = exp(bessel_k_scaled_log(x,lambda+1.0_dp)- &
            bessel_k_scaled_log(x,lambda))/x
      end if
   end function kappa_gh

   function delta_kappa_gh(x, lambda) result(value)
      real(dp), intent(in) :: x, lambda
      real(dp) :: value
      value = kappa_gh(x,lambda+1.0_dp)-kappa_gh(x,lambda)
   end function delta_kappa_gh

   subroutine gh_parameters(skew, shape, lambda, alpha, beta, delta, mu, valid)
      real(dp), intent(in) :: skew, shape, lambda
      real(dp), intent(out) :: alpha, beta, delta, mu
      logical, intent(out), optional :: valid
      real(dp) :: rho, rho2, kap, dkap, a2
      logical :: ok

      rho = max(-0.999999_dp,min(0.999999_dp,skew))
      rho2 = 1.0_dp-rho*rho
      ok = shape > 0.0_dp .and. rho2 > 0.0_dp
      if (.not. ok) then
         alpha=0.0_dp; beta=0.0_dp; delta=0.0_dp; mu=0.0_dp
         if (present(valid)) valid=.false.
         return
      end if
      kap = kappa_gh(shape,lambda)
      dkap = delta_kappa_gh(shape,lambda)
      a2 = shape*shape*kap/rho2
      a2 = a2*(1.0_dp+rho*rho*shape*shape*dkap/rho2)
      ok = a2 > 0.0_dp
      if (ok) then
         alpha = sqrt(a2)
         beta = alpha*rho
         delta = shape/(alpha*sqrt(rho2))
         mu = -beta*delta*delta*kap
      else
         alpha=0.0_dp; beta=0.0_dp; delta=0.0_dp; mu=0.0_dp
      end if
      if (present(valid)) valid=ok
   end subroutine gh_parameters

   subroutine ghst_parameters(skew, shape, beta, delta, mu, valid)
      real(dp), intent(in) :: skew, shape
      real(dp), intent(out) :: beta, delta, mu
      logical, intent(out), optional :: valid
      real(dp) :: den
      logical :: ok

      ok = shape > 4.0_dp
      if (ok) then
         den = 2.0_dp*skew*skew/((shape-2.0_dp)**2*(shape-4.0_dp)) + 1.0_dp/(shape-2.0_dp)
         ok = den > 0.0_dp
      end if
      if (ok) then
         delta = 1.0_dp/sqrt(den)
         beta = skew/delta
         mu = -beta*delta*delta/(shape-2.0_dp)
      else
         beta=0.0_dp; delta=0.0_dp; mu=0.0_dp
      end if
      if (present(valid)) valid=ok
   end subroutine ghst_parameters

   function dgh_raw(x, alpha, beta, delta, mu, lambda, log_density) result(value)
      real(dp), intent(in) :: x, alpha, beta, delta, mu, lambda
      logical, intent(in), optional :: log_density
      real(dp) :: value, arg0, argx, xm, logpdf
      logical :: want_log

      want_log=.false.
      if (present(log_density)) want_log=log_density
      if (alpha<=0.0_dp .or. delta<=0.0_dp .or. abs(beta)>=alpha) then
         value=merge(-huge(1.0_dp),0.0_dp,want_log)
         return
      end if
      xm=x-mu
      arg0=delta*sqrt(alpha*alpha-beta*beta)
      argx=alpha*sqrt(delta*delta+xm*xm)
      logpdf = 0.5_dp*lambda*log(alpha*alpha-beta*beta) - 0.5_dp*log(2.0_dp*pi) - &
         (lambda-0.5_dp)*log(alpha) - lambda*log(delta) - bessel_k_log(arg0,lambda) + &
         0.5_dp*(lambda-0.5_dp)*log(delta*delta+xm*xm) + &
         bessel_k_log(argx,lambda-0.5_dp) + beta*xm
      value=merge(logpdf,exp(max(log(tiny(1.0_dp)),min(log(huge(1.0_dp)),logpdf))),want_log)
   end function dgh_raw

   function dsgh(x, mean, sd, skew, shape, lambda, log_density) result(value)
      real(dp), intent(in) :: x, mean, sd, skew, shape, lambda
      logical, intent(in), optional :: log_density
      real(dp) :: value, alpha, beta, delta, mu, lp
      logical :: ok, want_log
      want_log=.false.; if(present(log_density)) want_log=log_density
      call gh_parameters(skew,shape,lambda,alpha,beta,delta,mu,ok)
      if (.not.ok .or. sd<=0.0_dp) then
         value=merge(-huge(1.0_dp),0.0_dp,want_log)
         return
      end if
      lp=dgh_raw((x-mean)/sd,alpha,beta,delta,mu,lambda,.true.)-log(sd)
      value=merge(lp,exp(max(log(tiny(1.0_dp)),min(log(huge(1.0_dp)),lp))),want_log)
   end function dsgh

   function dsnig(x, mean, sd, skew, shape, log_density) result(value)
      real(dp), intent(in) :: x, mean, sd, skew, shape
      logical, intent(in), optional :: log_density
      real(dp) :: value
      if (present(log_density)) then
         value=dsgh(x,mean,sd,skew,shape,-0.5_dp,log_density)
      else
         value=dsgh(x,mean,sd,skew,shape,-0.5_dp)
      end if
   end function dsnig

   function dsghst(x, mean, sd, skew, shape, log_density) result(value)
      real(dp), intent(in) :: x, mean, sd, skew, shape
      logical, intent(in), optional :: log_density
      real(dp) :: value, beta, delta, mu, z, xm, arg, lp, b
      logical :: ok, want_log
      want_log=.false.; if(present(log_density)) want_log=log_density
      if (sd<=0.0_dp) then
         value=merge(-huge(1.0_dp),0.0_dp,want_log); return
      end if
      call ghst_parameters(skew,shape,beta,delta,mu,ok)
      if (.not.ok) then
         value=merge(-huge(1.0_dp),0.0_dp,want_log); return
      end if
      z=(x-mean)/sd
      if (abs(beta)<1.0e-8_dp) then
         lp=log(student_t_pdf(z*sqrt(shape/(shape-2.0_dp)),shape)) + &
            0.5_dp*log(shape/(shape-2.0_dp))-log(sd)
      else
         xm=z-mu
         b=abs(beta)
         arg=b*sqrt(delta*delta+xm*xm)
         lp=0.5_dp*(1.0_dp-shape)*log(2.0_dp)+shape*log(delta)+ &
            0.5_dp*(shape+1.0_dp)*log(b)+bessel_k_log(arg,0.5_dp*(shape+1.0_dp))+ &
            beta*xm-log_gamma(0.5_dp*shape)-0.5_dp*log(pi)- &
            0.25_dp*(shape+1.0_dp)*log(delta*delta+xm*xm)-log(sd)
      end if
      value=merge(lp,exp(max(log(tiny(1.0_dp)),min(log(huge(1.0_dp)),lp))),want_log)
   end function dsghst

   function transformed_integrand(u, kind, skew, shape, lambda) result(value)
      real(dp), intent(in) :: u, skew, shape, lambda
      integer, intent(in) :: kind
      real(dp) :: value, x
      if (u<=1.0e-12_dp .or. u>=1.0_dp-1.0e-12_dp) then
         value=0.0_dp
         return
      end if
      x=tan(pi*(u-0.5_dp))
      select case(kind)
      case(1); value=dsgh(x,0.0_dp,1.0_dp,skew,shape,lambda)
      case(2); value=dsnig(x,0.0_dp,1.0_dp,skew,shape)
      case default; value=dsghst(x,0.0_dp,1.0_dp,skew,shape)
      end select
      value=value*pi*(1.0_dp+x*x)
   end function transformed_integrand

   function gh_cdf_numeric(q, kind, skew, shape, lambda) result(value)
      real(dp), intent(in) :: q, skew, shape, lambda
      integer, intent(in) :: kind
      real(dp) :: value, uq, h, u
      integer :: i, n
      uq=0.5_dp+atan(q)/pi
      if (uq<=1.0e-12_dp) then
         value=0.0_dp; return
      else if (uq>=1.0_dp-1.0e-12_dp) then
         value=1.0_dp; return
      end if
      n=max(160,2*ceiling(320.0_dp*uq))
      if(mod(n,2)/=0)n=n+1
      h=uq/real(n,dp)
      value=transformed_integrand(0.0_dp,kind,skew,shape,lambda)+ &
         transformed_integrand(uq,kind,skew,shape,lambda)
      do i=1,n-1
         u=real(i,dp)*h
         if(mod(i,2)==0) then
            value=value+2.0_dp*transformed_integrand(u,kind,skew,shape,lambda)
         else
            value=value+4.0_dp*transformed_integrand(u,kind,skew,shape,lambda)
         end if
      end do
      value=max(0.0_dp,min(1.0_dp,value*h/3.0_dp))
   end function gh_cdf_numeric

   function psgh(q, mean, sd, skew, shape, lambda) result(value)
      real(dp), intent(in) :: q, mean, sd, skew, shape, lambda
      real(dp) :: value
      if(sd<=0.0_dp) then; value=0.0_dp; else; value=gh_cdf_numeric((q-mean)/sd,1,skew,shape,lambda); end if
   end function psgh

   function psnig(q, mean, sd, skew, shape) result(value)
      real(dp), intent(in) :: q, mean, sd, skew, shape
      real(dp) :: value
      if(sd<=0.0_dp) then; value=0.0_dp; else; value=gh_cdf_numeric((q-mean)/sd,2,skew,shape,-0.5_dp); end if
   end function psnig

   function psghst(q, mean, sd, skew, shape) result(value)
      real(dp), intent(in) :: q, mean, sd, skew, shape
      real(dp) :: value
      if(sd<=0.0_dp) then; value=0.0_dp; else; value=gh_cdf_numeric((q-mean)/sd,3,skew,shape,-0.5_dp*shape); end if
   end function psghst

   function gh_quantile_numeric(p, kind, skew, shape, lambda) result(value)
      real(dp), intent(in) :: p, skew, shape, lambda
      integer, intent(in) :: kind
      real(dp) :: value, lo, hi, mid, c
      integer :: i
      if(p<=0.0_dp) then; value=-huge(1.0_dp); return; end if
      if(p>=1.0_dp) then; value=huge(1.0_dp); return; end if
      lo=-2.0_dp; hi=2.0_dp
      do while(gh_cdf_numeric(lo,kind,skew,shape,lambda)>p .and. abs(lo)<1.0e8_dp); lo=2.0_dp*lo; end do
      do while(gh_cdf_numeric(hi,kind,skew,shape,lambda)<p .and. abs(hi)<1.0e8_dp); hi=2.0_dp*hi; end do
      do i=1,44
         mid=0.5_dp*(lo+hi)
         c=gh_cdf_numeric(mid,kind,skew,shape,lambda)
         if(c<p)then; lo=mid; else; hi=mid; end if
      end do
      value=0.5_dp*(lo+hi)
   end function gh_quantile_numeric

   function qsgh(p, mean, sd, skew, shape, lambda) result(value)
      real(dp), intent(in) :: p, mean, sd, skew, shape, lambda
      real(dp) :: value
      value=mean+sd*gh_quantile_numeric(p,1,skew,shape,lambda)
   end function qsgh

   function qsnig(p, mean, sd, skew, shape) result(value)
      real(dp), intent(in) :: p, mean, sd, skew, shape
      real(dp) :: value
      value=mean+sd*gh_quantile_numeric(p,2,skew,shape,-0.5_dp)
   end function qsnig

   function qsghst(p, mean, sd, skew, shape) result(value)
      real(dp), intent(in) :: p, mean, sd, skew, shape
      real(dp) :: value
      value=mean+sd*gh_quantile_numeric(p,3,skew,shape,-0.5_dp*shape)
   end function qsghst

   pure function gig_g(y,m,beta,lambda) result(value)
      real(dp), intent(in)::y,m,beta,lambda
      real(dp)::value
      value=0.5_dp*beta*y**3-y*y*(0.5_dp*beta*m+lambda+1.0_dp)+ &
         y*((lambda-1.0_dp)*m-0.5_dp*beta)+0.5_dp*beta*m
   end function gig_g

   function gig_root(lo0,hi0,m,beta,lambda) result(root)
      real(dp),intent(in)::lo0,hi0,m,beta,lambda
      real(dp)::root,lo,hi,mid,flo,fmid
      integer::i
      lo=max(lo0,1.0e-14_dp);hi=hi0;flo=gig_g(lo,m,beta,lambda)
      do i=1,100
         mid=0.5_dp*(lo+hi);fmid=gig_g(mid,m,beta,lambda)
         if(flo*fmid<=0.0_dp)then;hi=mid;else;lo=mid;flo=fmid;end if
      end do
      root=0.5_dp*(lo+hi)
   end function gig_root

   function random_gig(lambda, chi, psi) result(value)
      real(dp), intent(in) :: lambda,chi,psi
      real(dp) :: value,alpha,beta,m,upper,ym,yp,a,b,c,r1,r2,y,lm1,m1,u
      integer :: iter
      if(chi<1.0e-14_dp .and. psi<1.0e-14_dp)then;value=-1.0_dp;return;end if
      if(chi<1.0e-14_dp)then
         if(lambda>0.0_dp)then;value=random_gamma(lambda)*2.0_dp/psi;else;value=-1.0_dp;end if
         return
      end if
      if(psi<1.0e-14_dp)then
         if(lambda<0.0_dp)then;value=1.0_dp/(random_gamma(-lambda)*2.0_dp/chi);else;value=-1.0_dp;end if
         return
      end if
      alpha=sqrt(psi/chi);beta=sqrt(psi*chi);lm1=lambda-1.0_dp
      m=(lm1+sqrt(lm1*lm1+beta*beta))/beta
      m1=m+1.0_dp/m;upper=m
      do while(gig_g(upper,m,beta,lambda)<=0.0_dp .and. upper<1.0e12_dp);upper=2.0_dp*upper;end do
      ym=gig_root(0.0_dp,m,m,beta,lambda);yp=gig_root(m,upper,m,beta,lambda)
      a=(yp-m)*(yp/m)**(0.5_dp*lm1)*exp(-0.25_dp*beta*(yp+1.0_dp/yp-m1))
      b=(ym-m)*(ym/m)**(0.5_dp*lm1)*exp(-0.25_dp*beta*(ym+1.0_dp/ym-m1))
      c=-0.25_dp*beta*m1+0.5_dp*lm1*log(m)
      y=m
      do iter=1,100000
         call random_number(r1);call random_number(r2);r1=max(r1,tiny(1.0_dp))
         y=m+a*r2/r1+b*(1.0_dp-r2)/r1
         if(y>0.0_dp)then
            if(-log(r1)>=-0.5_dp*lm1*log(y)+0.25_dp*beta*(y+1.0_dp/y)+c)exit
         end if
      end do
      if(iter>100000)then
         call random_number(u); y=max(m*(0.5_dp+u),tiny(1.0_dp))
      end if
      value=y/alpha
   end function random_gig

   function rsgh(mean,sd,skew,shape,lambda) result(value)
      real(dp),intent(in)::mean,sd,skew,shape,lambda
      real(dp)::value,alpha,beta,delta,mu,w
      logical::ok
      call gh_parameters(skew,shape,lambda,alpha,beta,delta,mu,ok)
      if(.not.ok)then;value=0.0_dp;return;end if
      w=random_gig(lambda,delta*delta,alpha*alpha-beta*beta)
      value=mean+sd*(mu+beta*w+sqrt(max(w,0.0_dp))*random_normal())
   end function rsgh

   function rsnig(mean,sd,skew,shape) result(value)
      real(dp),intent(in)::mean,sd,skew,shape
      real(dp)::value
      value=rsgh(mean,sd,skew,shape,-0.5_dp)
   end function rsnig

   function rsghst(mean,sd,skew,shape) result(value)
      real(dp),intent(in)::mean,sd,skew,shape
      real(dp)::value,beta,delta,mu,w
      logical::ok
      call ghst_parameters(skew,shape,beta,delta,mu,ok)
      if(.not.ok)then;value=0.0_dp;return;end if
      w=1.0_dp/(random_gamma(0.5_dp*shape)*2.0_dp/(delta*delta))
      value=mean+sd*(mu+beta*w+sqrt(w)*random_normal())
   end function rsghst

end module rugarch_gh
