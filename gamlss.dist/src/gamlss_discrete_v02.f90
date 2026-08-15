! Extended discrete families from gamlss.dist 6.1-1.
! Original package GPL-2 | GPL-3. Translation SPDX-License-Identifier: GPL-3.0-only
module gamlss_discrete_v02
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf
   use gamlss_kinds, only : dp
   use gamlss_special, only : log_beta_fn, zeta_fn
   use gamlss_base, only : dpois_v, ppois_v
   use gamlss_v02_numerics, only : log_bessel_k
   implicit none
   private
   public :: dGPO,pGPO,qGPO,rGPO,dDPO,pDPO,qDPO,rDPO
   public :: dDEL,pDEL,qDEL,rDEL,dSI,pSI,qSI,rSI,dSICHEL,pSICHEL,qSICHEL,rSICHEL
   public :: dYULE,pYULE,qYULE,rYULE,dWARING,pWARING,qWARING,rWARING
   public :: dZIPF,pZIPF,qZIPF,rZIPF
contains

   elemental real(dp) function nanv() result(x)
      x=ieee_value(0.0_dp,ieee_quiet_nan)
   end function nanv

   elemental real(dp) function infv() result(x)
      x=ieee_value(0.0_dp,ieee_positive_inf)
   end function infv

   elemental logical function want_log(flag) result(v)
      logical, intent(in), optional :: flag
      v=.false.
      if (present(flag)) v=flag
   end function want_log

   elemental logical function lower(flag) result(v)
      logical, intent(in), optional :: flag
      v=.true.
      if (present(flag)) v=flag
   end function lower

   elemental logical function is_int(x) result(ok)
      real(dp), intent(in) :: x
      ok=abs(x-real(nint(x),dp))<1.0e-10_dp
   end function is_int

   elemental real(dp) function finish_prob(p,lower_tail,log_p) result(v)
      real(dp), intent(in) :: p
      logical, intent(in), optional :: lower_tail,log_p
      v=max(0.0_dp,min(1.0_dp,p))
      if (.not.lower(lower_tail)) v=1.0_dp-v
      if (present(log_p)) then
         if (log_p) then
            if (v<=0.0_dp) then
               v=-infv()
            else
               v=log(v)
            end if
         end if
      end if
   end function finish_prob

   elemental real(dp) function input_prob(p,lower_tail,log_p) result(v)
      real(dp), intent(in) :: p
      logical, intent(in), optional :: lower_tail,log_p
      v=p
      if (present(log_p)) then
         if (log_p) v=exp(v)
      end if
      if (.not.lower(lower_tail)) v=1.0_dp-v
   end function input_prob

   elemental real(dp) function dGPO(x,mu,sigma,log_density) result(v)
      real(dp), intent(in) :: x,mu,sigma
      logical, intent(in), optional :: log_density
      real(dp) :: ld
      integer :: k
      if (mu<=0.0_dp .or. sigma<0.0_dp) then
         v=nanv()
      else if (x<0.0_dp .or. .not.is_int(x)) then
         v=merge(-infv(),0.0_dp,want_log(log_density))
      else if (sigma<1.0e-8_dp) then
         v=dpois_v(nint(x),mu,want_log(log_density))
      else
         k=nint(x)
         ld=real(k,dp)*log(mu/(1.0_dp+sigma*mu))
         if (k>0) ld=ld+real(k-1,dp)*log(1.0_dp+sigma*real(k,dp))
         ld=ld-mu*(1.0_dp+sigma*real(k,dp))/(1.0_dp+sigma*mu)-log_gamma(real(k+1,dp))
         v=merge(ld,exp(ld),want_log(log_density))
      end if
   end function dGPO

   real(dp) function pGPO(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,mu,sigma
      logical, intent(in), optional :: lower_tail,log_p
      integer :: k,j
      real(dp) :: s
      if (mu<=0.0_dp .or. sigma<0.0_dp) then
         v=nanv()
         return
      end if
      k=int(floor(q))
      if (k<0) then
         s=0.0_dp
      else if (sigma<1.0e-8_dp) then
         s=ppois_v(k,mu)
      else
         s=0.0_dp
         do j=0,k
            s=s+dGPO(real(j,dp),mu,sigma)
         end do
      end if
      v=finish_prob(s,lower_tail,log_p)
   end function pGPO

   integer function qGPO(p,mu,sigma,lower_tail,log_p,max_value) result(v)
      real(dp), intent(in) :: p,mu,sigma
      logical, intent(in), optional :: lower_tail,log_p
      integer, intent(in), optional :: max_value
      real(dp) :: pp,cum
      integer :: j,mx
      pp=input_prob(p,lower_tail,log_p)
      mx=10000
      if (present(max_value)) mx=max_value
      if (pp<0.0_dp .or. pp>1.0_dp .or. mu<=0.0_dp .or. sigma<0.0_dp) then
         v=-huge(v)
      else if (pp>=1.0_dp) then
         v=huge(v)
      else
         cum=0.0_dp
         v=mx
         do j=0,mx
            cum=cum+dGPO(real(j,dp),mu,sigma)
            if (cum>=pp) then
               v=j
               exit
            end if
         end do
      end if
   end function qGPO

   integer function rGPO(mu,sigma,max_value) result(v)
      real(dp), intent(in) :: mu,sigma
      integer, intent(in), optional :: max_value
      real(dp) :: u
      call random_number(u)
      v=qGPO(u,mu,sigma,max_value=max_value)
   end function rGPO

   real(dp) function dpo_log_kernel(k,mu,sigma) result(lh)
      integer, intent(in) :: k
      real(dp), intent(in) :: mu,sigma
      real(dp) :: y,yl
      y=real(k,dp)
      yl=0.0_dp
      if (k>0) yl=log(y)
      lh=-0.5_dp*log(sigma)-mu/sigma-log_gamma(y+1.0_dp)
      lh=lh+y*yl-y+(y*log(mu))/sigma+y/sigma-(y*yl)/sigma
   end function dpo_log_kernel

   real(dp) function dpo_log_c(mu,sigma,maxk) result(logc)
      real(dp), intent(in) :: mu,sigma
      integer, intent(in) :: maxk
      real(dp) :: m,s,lk
      integer :: k
      m=-huge(1.0_dp)
      do k=0,maxk
         lk=dpo_log_kernel(k,mu,sigma)
         m=max(m,lk)
      end do
      s=0.0_dp
      do k=0,maxk
         s=s+exp(dpo_log_kernel(k,mu,sigma)-m)
      end do
      logc=-m-log(s)
   end function dpo_log_c

   real(dp) function dDPO(x,mu,sigma,log_density) result(v)
      real(dp), intent(in) :: x,mu,sigma
      logical, intent(in), optional :: log_density
      integer :: k,maxk
      real(dp) :: ld
      if (mu<=0.0_dp .or. sigma<=0.0_dp) then
         v=nanv()
      else if (x<0.0_dp .or. .not.is_int(x)) then
         v=merge(-infv(),0.0_dp,want_log(log_density))
      else
         k=nint(x)
         maxk=max(500,3*k)
         ld=dpo_log_kernel(k,mu,sigma)+dpo_log_c(mu,sigma,maxk)
         v=merge(ld,exp(ld),want_log(log_density))
      end if
   end function dDPO

   real(dp) function pDPO(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,mu,sigma
      logical, intent(in), optional :: lower_tail,log_p
      integer :: k,j,maxk
      real(dp) :: s,logc
      if (mu<=0.0_dp .or. sigma<=0.0_dp) then
         v=nanv()
         return
      end if
      k=int(floor(q))
      if (k<0) then
         s=0.0_dp
      else
         maxk=max(500,3*k)
         logc=dpo_log_c(mu,sigma,maxk)
         s=0.0_dp
         do j=0,k
            s=s+exp(dpo_log_kernel(j,mu,sigma)+logc)
         end do
      end if
      v=finish_prob(s,lower_tail,log_p)
   end function pDPO

   integer function qDPO(p,mu,sigma,lower_tail,log_p,max_value) result(v)
      real(dp), intent(in) :: p,mu,sigma
      logical, intent(in), optional :: lower_tail,log_p
      integer, intent(in), optional :: max_value
      real(dp) :: pp,cum,logc
      integer :: j,mx,maxk
      pp=input_prob(p,lower_tail,log_p)
      mx=10000
      if (present(max_value)) mx=max_value
      if (pp<0.0_dp .or. pp>1.0_dp .or. mu<=0.0_dp .or. sigma<=0.0_dp) then
         v=-huge(v)
      else if (pp>=1.0_dp) then
         v=huge(v)
      else
         maxk=max(500,3*mx)
         logc=dpo_log_c(mu,sigma,maxk)
         cum=0.0_dp
         v=mx
         do j=0,mx
            cum=cum+exp(dpo_log_kernel(j,mu,sigma)+logc)
            if (cum>=pp) then
               v=j
               exit
            end if
         end do
      end if
   end function qDPO

   integer function rDPO(mu,sigma,max_value) result(v)
      real(dp), intent(in) :: mu,sigma
      integer, intent(in), optional :: max_value
      real(dp) :: u
      call random_number(u)
      v=qDPO(u,mu,sigma,max_value=max_value)
   end function rDPO

   real(dp) function dDEL(x,mu,sigma,nu,log_density) result(v)
      real(dp), intent(in) :: x,mu,sigma,nu
      logical, intent(in), optional :: log_density
      integer :: k,j
      real(dp) :: logp0,tprev,tnew,sumlog,dum,ld
      if (mu<=0.0_dp .or. sigma<0.0_dp .or. nu<=0.0_dp .or. nu>=1.0_dp) then
         v=nanv()
      else if (x<0.0_dp .or. .not.is_int(x)) then
         v=merge(-infv(),0.0_dp,want_log(log_density))
      else if (sigma<1.0e-8_dp) then
         v=dpois_v(nint(x),mu,want_log(log_density))
      else
         k=nint(x)
         logp0=-mu*nu-log(1.0_dp+mu*sigma*(1.0_dp-nu))/sigma
         tprev=mu*nu+mu*(1.0_dp-nu)/(1.0_dp+mu*sigma*(1.0_dp-nu))
         sumlog=0.0_dp
         do j=1,k
            if (tprev<=0.0_dp) then
               v=nanv()
               return
            end if
            sumlog=sumlog+log(tprev)
            dum=1.0_dp+1.0_dp/(mu*sigma*(1.0_dp-nu))
            tnew=(real(j,dp)+mu*nu+1.0_dp/(sigma*(1.0_dp-nu)) &
               -mu*nu*real(j,dp)/tprev)/dum
            tprev=tnew
         end do
         ld=logp0-log_gamma(real(k+1,dp))+sumlog
         v=merge(ld,exp(ld),want_log(log_density))
      end if
   end function dDEL

   real(dp) function pDEL(q,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,mu,sigma,nu
      logical, intent(in), optional :: lower_tail,log_p
      integer :: k,j
      real(dp) :: s
      if (mu<=0.0_dp .or. sigma<0.0_dp .or. nu<=0.0_dp .or. nu>=1.0_dp) then
         v=nanv()
         return
      end if
      k=int(floor(q))
      s=0.0_dp
      if (k>=0) then
         do j=0,k
            s=s+dDEL(real(j,dp),mu,sigma,nu)
         end do
      end if
      v=finish_prob(s,lower_tail,log_p)
   end function pDEL

   integer function qDEL(p,mu,sigma,nu,lower_tail,log_p,max_value) result(v)
      real(dp), intent(in) :: p,mu,sigma,nu
      logical, intent(in), optional :: lower_tail,log_p
      integer, intent(in), optional :: max_value
      real(dp) :: pp,cum
      integer :: j,mx
      pp=input_prob(p,lower_tail,log_p)
      mx=10000
      if (present(max_value)) mx=max_value
      if (pp<0.0_dp .or. pp>1.0_dp) then
         v=-huge(v)
      else if (pp>=1.0_dp) then
         v=huge(v)
      else
         cum=0.0_dp
         v=mx
         do j=0,mx
            cum=cum+dDEL(real(j,dp),mu,sigma,nu)
            if (cum>=pp) then
               v=j
               exit
            end if
         end do
      end if
   end function qDEL

   integer function rDEL(mu,sigma,nu,max_value) result(v)
      real(dp), intent(in) :: mu,sigma,nu
      integer, intent(in), optional :: max_value
      real(dp) :: u
      call random_number(u)
      v=qDEL(u,mu,sigma,nu,max_value=max_value)
   end function rDEL

   real(dp) function dSI(x,mu,sigma,nu,log_density) result(v)
      real(dp), intent(in) :: x,mu,sigma,nu
      logical, intent(in), optional :: log_density
      integer :: k,j
      real(dp) :: alpha,lbes,tprev,tnew,sumlog,ld
      if (mu<=0.0_dp .or. sigma<=0.0_dp) then
         v=nanv()
      else if (x<0.0_dp .or. .not.is_int(x)) then
         v=merge(-infv(),0.0_dp,want_log(log_density))
      else
         k=nint(x)
         alpha=sqrt(1.0_dp+2.0_dp*sigma*mu)/sigma
         lbes=log_bessel_k(alpha,nu+1.0_dp)-log_bessel_k(alpha,nu)
         tprev=mu/sqrt(1.0_dp+2.0_dp*sigma*mu)*exp(lbes)
         sumlog=0.0_dp
         do j=1,k
            if (tprev<=0.0_dp) then
               v=nanv()
               return
            end if
            sumlog=sumlog+log(tprev)
            tnew=(sigma*2.0_dp*(real(j,dp)+nu)/mu+1.0_dp/tprev) &
               *(mu/(sigma*alpha))**2
            tprev=tnew
         end do
         ld=-log_gamma(real(k+1,dp))-nu*log(sigma*alpha)+sumlog
         ld=ld+log_bessel_k(alpha,nu)-log_bessel_k(1.0_dp/sigma,nu)
         v=merge(ld,exp(ld),want_log(log_density))
      end if
   end function dSI

   real(dp) function pSI(q,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,mu,sigma,nu
      logical, intent(in), optional :: lower_tail,log_p
      integer :: k,j
      real(dp) :: s
      if (mu<=0.0_dp .or. sigma<=0.0_dp) then
         v=nanv()
         return
      end if
      k=int(floor(q))
      s=0.0_dp
      if (k>=0) then
         do j=0,k
            s=s+dSI(real(j,dp),mu,sigma,nu)
         end do
      end if
      v=finish_prob(s,lower_tail,log_p)
   end function pSI

   integer function qSI(p,mu,sigma,nu,lower_tail,log_p,max_value) result(v)
      real(dp), intent(in) :: p,mu,sigma,nu
      logical, intent(in), optional :: lower_tail,log_p
      integer, intent(in), optional :: max_value
      real(dp) :: pp,cum
      integer :: j,mx
      pp=input_prob(p,lower_tail,log_p)
      mx=10000
      if (present(max_value)) mx=max_value
      if (pp<0.0_dp .or. pp>1.0_dp) then
         v=-huge(v)
      else if (pp>=1.0_dp) then
         v=huge(v)
      else
         cum=0.0_dp
         v=mx
         do j=0,mx
            cum=cum+dSI(real(j,dp),mu,sigma,nu)
            if (cum>=pp) then
               v=j
               exit
            end if
         end do
      end if
   end function qSI

   integer function rSI(mu,sigma,nu,max_value) result(v)
      real(dp), intent(in) :: mu,sigma,nu
      integer, intent(in), optional :: max_value
      real(dp) :: u
      call random_number(u)
      v=qSI(u,mu,sigma,nu,max_value=max_value)
   end function rSI

   real(dp) function dSICHEL(x,mu,sigma,nu,log_density) result(v)
      real(dp), intent(in) :: x,mu,sigma,nu
      logical, intent(in), optional :: log_density
      integer :: k,j
      real(dp) :: cvec,alpha,lbes,tprev,tnew,sumlog,ld
      if (mu<=0.0_dp .or. sigma<=0.0_dp) then
         v=nanv()
      else if (x<0.0_dp .or. .not.is_int(x)) then
         v=merge(-infv(),0.0_dp,want_log(log_density))
      else
         k=nint(x)
         cvec=exp(log_bessel_k(1.0_dp/sigma,nu+1.0_dp)-log_bessel_k(1.0_dp/sigma,nu))
         alpha=sqrt(1.0_dp+2.0_dp*sigma*mu/cvec)/sigma
         lbes=log_bessel_k(alpha,nu+1.0_dp)-log_bessel_k(alpha,nu)
         tprev=(mu/cvec)/sqrt(1.0_dp+2.0_dp*sigma*mu/cvec)*exp(lbes)
         sumlog=0.0_dp
         do j=1,k
            if (tprev<=0.0_dp) then
               v=nanv()
               return
            end if
            sumlog=sumlog+log(tprev)
            tnew=(cvec*sigma*2.0_dp*(real(j,dp)+nu)/mu+1.0_dp/tprev) &
               *(mu/(sigma*alpha*cvec))**2
            tprev=tnew
         end do
         ld=-log_gamma(real(k+1,dp))-nu*log(sigma*alpha)+sumlog
         ld=ld+log_bessel_k(alpha,nu)-log_bessel_k(1.0_dp/sigma,nu)
         v=merge(ld,exp(ld),want_log(log_density))
      end if
   end function dSICHEL

   real(dp) function pSICHEL(q,mu,sigma,nu,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,mu,sigma,nu
      logical, intent(in), optional :: lower_tail,log_p
      integer :: k,j
      real(dp) :: s
      if (mu<=0.0_dp .or. sigma<=0.0_dp) then
         v=nanv()
         return
      end if
      k=int(floor(q))
      s=0.0_dp
      if (k>=0) then
         do j=0,k
            s=s+dSICHEL(real(j,dp),mu,sigma,nu)
         end do
      end if
      v=finish_prob(s,lower_tail,log_p)
   end function pSICHEL

   integer function qSICHEL(p,mu,sigma,nu,lower_tail,log_p,max_value) result(v)
      real(dp), intent(in) :: p,mu,sigma,nu
      logical, intent(in), optional :: lower_tail,log_p
      integer, intent(in), optional :: max_value
      real(dp) :: pp,cum
      integer :: j,mx
      pp=input_prob(p,lower_tail,log_p)
      mx=10000
      if (present(max_value)) mx=max_value
      if (pp<0.0_dp .or. pp>1.0_dp) then
         v=-huge(v)
      else if (pp>=1.0_dp) then
         v=huge(v)
      else
         cum=0.0_dp
         v=mx
         do j=0,mx
            cum=cum+dSICHEL(real(j,dp),mu,sigma,nu)
            if (cum>=pp) then
               v=j
               exit
            end if
         end do
      end if
   end function qSICHEL

   integer function rSICHEL(mu,sigma,nu,max_value) result(v)
      real(dp), intent(in) :: mu,sigma,nu
      integer, intent(in), optional :: max_value
      real(dp) :: u
      call random_number(u)
      v=qSICHEL(u,mu,sigma,nu,max_value=max_value)
   end function rSICHEL

   elemental real(dp) function dYULE(x,mu,log_density) result(v)
      real(dp), intent(in) :: x,mu
      logical, intent(in), optional :: log_density
      real(dp) :: lambda,ld
      if (mu<=0.0_dp) then
         v=nanv()
      else if (x<0.0_dp .or. .not.is_int(x)) then
         v=merge(-infv(),0.0_dp,want_log(log_density))
      else
         lambda=(mu+1.0_dp)/mu
         ld=log_beta_fn(lambda+1.0_dp,x+1.0_dp)-log_beta_fn(lambda,1.0_dp)
         v=merge(ld,exp(ld),want_log(log_density))
      end if
   end function dYULE

   real(dp) function pYULE(q,mu,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,mu
      logical, intent(in), optional :: lower_tail,log_p
      integer :: k,j
      real(dp) :: s
      k=int(floor(q))
      s=0.0_dp
      if (k>=0) then
         do j=0,k
            s=s+dYULE(real(j,dp),mu)
         end do
      end if
      v=finish_prob(s,lower_tail,log_p)
   end function pYULE

   integer function qYULE(p,mu,lower_tail,log_p,max_value) result(v)
      real(dp), intent(in) :: p,mu
      logical, intent(in), optional :: lower_tail,log_p
      integer, intent(in), optional :: max_value
      real(dp) :: pp,cum
      integer :: j,mx
      pp=input_prob(p,lower_tail,log_p)
      mx=10000
      if (present(max_value)) mx=max_value
      cum=0.0_dp
      v=mx
      if (pp>=1.0_dp) then
         v=huge(v)
         return
      end if
      do j=0,mx
         cum=cum+dYULE(real(j,dp),mu)
         if (cum>=pp) then
            v=j
            exit
         end if
      end do
   end function qYULE

   integer function rYULE(mu,max_value) result(v)
      real(dp), intent(in) :: mu
      integer, intent(in), optional :: max_value
      real(dp) :: u
      call random_number(u)
      v=qYULE(u,mu,max_value=max_value)
   end function rYULE

   elemental real(dp) function dWARING(x,mu,sigma,log_density) result(v)
      real(dp), intent(in) :: x,mu,sigma
      logical, intent(in), optional :: log_density
      real(dp) :: ld
      if (mu<=0.0_dp .or. sigma<=0.0_dp) then
         v=nanv()
      else if (x<0.0_dp .or. .not.is_int(x)) then
         v=merge(-infv(),0.0_dp,want_log(log_density))
      else
         ld=log_beta_fn(x+mu/sigma,1.0_dp/sigma+2.0_dp)
         ld=ld-log_beta_fn(mu/sigma,1.0_dp/sigma+1.0_dp)
         v=merge(ld,exp(ld),want_log(log_density))
      end if
   end function dWARING

   real(dp) function pWARING(q,mu,sigma,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,mu,sigma
      logical, intent(in), optional :: lower_tail,log_p
      integer :: k,j
      real(dp) :: s
      k=int(floor(q))
      s=0.0_dp
      if (k>=0) then
         do j=0,k
            s=s+dWARING(real(j,dp),mu,sigma)
         end do
      end if
      v=finish_prob(s,lower_tail,log_p)
   end function pWARING

   integer function qWARING(p,mu,sigma,lower_tail,log_p,max_value) result(v)
      real(dp), intent(in) :: p,mu,sigma
      logical, intent(in), optional :: lower_tail,log_p
      integer, intent(in), optional :: max_value
      real(dp) :: pp,cum
      integer :: j,mx
      pp=input_prob(p,lower_tail,log_p)
      mx=10000
      if (present(max_value)) mx=max_value
      cum=0.0_dp
      v=mx
      if (pp>=1.0_dp) then
         v=huge(v)
         return
      end if
      do j=0,mx
         cum=cum+dWARING(real(j,dp),mu,sigma)
         if (cum>=pp) then
            v=j
            exit
         end if
      end do
   end function qWARING

   integer function rWARING(mu,sigma,max_value) result(v)
      real(dp), intent(in) :: mu,sigma
      integer, intent(in), optional :: max_value
      real(dp) :: u
      call random_number(u)
      v=qWARING(u,mu,sigma,max_value=max_value)
   end function rWARING

   elemental real(dp) function dZIPF(x,mu,log_density) result(v)
      real(dp), intent(in) :: x,mu
      logical, intent(in), optional :: log_density
      real(dp) :: ld,z
      if (mu<=0.0_dp) then
         v=nanv()
      else if (x<1.0_dp .or. .not.is_int(x)) then
         v=merge(-infv(),0.0_dp,want_log(log_density))
      else
         z=zeta_fn(mu+1.0_dp)
         ld=-(mu+1.0_dp)*log(x)-log(z)
         v=merge(ld,exp(ld),want_log(log_density))
      end if
   end function dZIPF

   real(dp) function pZIPF(q,mu,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,mu
      logical, intent(in), optional :: lower_tail,log_p
      integer :: k,j
      real(dp) :: s,z
      if (mu<=0.0_dp) then
         v=nanv()
         return
      end if
      k=int(floor(q))
      if (k<1) then
         s=0.0_dp
      else
         z=zeta_fn(mu+1.0_dp)
         s=0.0_dp
         do j=1,k
            s=s+real(j,dp)**(-mu-1.0_dp)/z
         end do
      end if
      v=finish_prob(s,lower_tail,log_p)
   end function pZIPF

   integer function qZIPF(p,mu,lower_tail,log_p,max_value) result(v)
      real(dp), intent(in) :: p,mu
      logical, intent(in), optional :: lower_tail,log_p
      integer, intent(in), optional :: max_value
      real(dp) :: pp,cum,z
      integer :: j,mx
      pp=input_prob(p,lower_tail,log_p)
      mx=10000
      if (present(max_value)) mx=max_value
      if (pp>=1.0_dp) then
         v=huge(v)
         return
      end if
      z=zeta_fn(mu+1.0_dp)
      cum=0.0_dp
      v=mx
      do j=1,mx
         cum=cum+real(j,dp)**(-mu-1.0_dp)/z
         if (cum>=pp) then
            v=j
            exit
         end if
      end do
   end function qZIPF

   integer function rZIPF(mu,max_value) result(v)
      real(dp), intent(in) :: mu
      integer, intent(in), optional :: max_value
      real(dp) :: u
      call random_number(u)
      v=qZIPF(u,mu,max_value=max_value)
   end function rZIPF

end module gamlss_discrete_v02
