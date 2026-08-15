! rmutil computational translation
! Copyright (C) 1998-2001 J.K. Lindsey
! Copyright (C) 2026 OpenAI (modern Fortran translation)
! SPDX-License-Identifier: GPL-2.0-or-later
module rmutil_discrete
   use rmutil_kinds, only : dp
   use rmutil_special, only : log_beta, log_choose, regularized_gamma_p, &
      regularized_gamma_q, negative_binomial_pmf, negative_binomial_cdf
   implicit none
   private
   public :: pbetabinom, dbetabinom, qbetabinom, rbetabinom
   public :: pdoublebinom, ddoublebinom, qdoublebinom, rdoublebinom
   public :: pmultbinom, dmultbinom, qmultbinom, rmultbinom
   public :: pdoublepois, ddoublepois, qdoublepois, rdoublepois
   public :: pmultpois, dmultpois, qmultpois, rmultpois
   public :: ppvfpois, dpvfpois, qpvfpois, rpvfpois
   public :: pgammacount, dgammacount, qgammacount, rgammacount
   public :: pconsul, dconsul, qconsul, rconsul
contains

   elemental real(dp) function dbetabinom(y,size,m,s,log_p) result(v)
      integer,intent(in)::y,size
      real(dp),intent(in)::m,s
      logical,intent(in),optional::log_p
      real(dp)::t,u,lv; logical::lg
      lg=.false.; if(present(log_p))lg=log_p
      if(y<0 .or. y>size)then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
      t=s*m; u=s*(1.0_dp-m)
      lv=log_beta(real(y,dp)+t,real(size-y,dp)+u)-log_beta(t,u)+log_choose(size,y)
      v=merge(lv,exp(lv),lg)
   end function dbetabinom

   real(dp) function pbetabinom(q,size,m,s) result(p)
      integer,intent(in)::q,size; real(dp),intent(in)::m,s
      integer::k
      if(q<0)then; p=0.0_dp; return; else if(q>=size)then; p=1.0_dp; return; end if
      p=0.0_dp; do k=0,q; p=p+dbetabinom(k,size,m,s); end do; p=min(1.0_dp,p)
   end function pbetabinom

   integer function qbetabinom(p,size,m,s) result(q)
      real(dp),intent(in)::p,m,s; integer,intent(in)::size
      q=discrete_quantile_bounded(p,size,cdf)
   contains
      function cdf(k) result(v); integer,intent(in)::k; real(dp)::v; v=pbetabinom(k,size,m,s); end function cdf
   end function qbetabinom

   function rbetabinom(n,size,m,s) result(x)
      integer,intent(in)::n,size; real(dp),intent(in)::m,s
      integer,allocatable::x(:); real(dp)::u; integer::i
      allocate(x(n)); do i=1,n; call random_number(u); x(i)=qbetabinom(u,size,m,s); end do
   end function rbetabinom

   real(dp) function double_binom_log_weight(y,n,m,s) result(lw)
      integer,intent(in)::y,n; real(dp),intent(in)::m,s
      lw=log_choose(n,y)+real(n,dp)*(s-1.0_dp)*log(real(max(n,1),dp)) + &
         real(y,dp)*s*log(m)+real(n-y,dp)*s*log(1.0_dp-m)
      if(y>0)lw=lw-real(y,dp)*(s-1.0_dp)*log(real(y,dp))
      if(y<n)lw=lw-real(n-y,dp)*(s-1.0_dp)*log(real(n-y,dp))
   end function double_binom_log_weight

   real(dp) function double_binom_log_norm(n,m,s) result(ln)
      integer,intent(in)::n; real(dp),intent(in)::m,s
      real(dp)::mx,sumv,lw; integer::k
      mx=-huge(1.0_dp)
      do k=0,n; mx=max(mx,double_binom_log_weight(k,n,m,s)); end do
      sumv=0.0_dp
      do k=0,n; lw=double_binom_log_weight(k,n,m,s); sumv=sumv+exp(lw-mx); end do
      ln=mx+log(sumv)
   end function double_binom_log_norm

   real(dp) function ddoublebinom(y,size,m,s,log_p) result(v)
      integer,intent(in)::y,size; real(dp),intent(in)::m,s
      logical,intent(in),optional::log_p
      logical::lg; real(dp)::lv
      lg=.false.; if(present(log_p))lg=log_p
      if(y<0 .or. y>size)then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
      lv=double_binom_log_weight(y,size,m,s)-double_binom_log_norm(size,m,s)
      v=merge(lv,exp(lv),lg)
   end function ddoublebinom

   real(dp) function pdoublebinom(q,size,m,s) result(p)
      integer,intent(in)::q,size; real(dp),intent(in)::m,s
      real(dp)::ln; integer::k,qq
      if(q<0)then; p=0.0_dp; return; end if
      qq=min(q,size); ln=double_binom_log_norm(size,m,s); p=0.0_dp
      do k=0,qq; p=p+exp(double_binom_log_weight(k,size,m,s)-ln); end do
      p=min(1.0_dp,p)
   end function pdoublebinom

   integer function qdoublebinom(p,size,m,s) result(q)
      real(dp),intent(in)::p,m,s; integer,intent(in)::size
      q=discrete_quantile_bounded(p,size,cdf)
   contains
      function cdf(k) result(v); integer,intent(in)::k; real(dp)::v; v=pdoublebinom(k,size,m,s); end function cdf
   end function qdoublebinom

   function rdoublebinom(n,size,m,s) result(x)
      integer,intent(in)::n,size; real(dp),intent(in)::m,s
      integer,allocatable::x(:); real(dp)::u; integer::i
      allocate(x(n)); do i=1,n; call random_number(u); x(i)=qdoublebinom(u,size,m,s); end do
   end function rdoublebinom

   real(dp) function mult_binom_log_weight_cdf(y,n,m,s) result(lw)
      integer,intent(in)::y,n; real(dp),intent(in)::m,s
      real(dp)::ss
      ss=log(s)
      lw=log_choose(n,y)+real(n-y,dp)*log(1.0_dp-m)+ &
         real(y,dp)*(log(m)+real((n-y)*y,dp)*ss)
   end function mult_binom_log_weight_cdf

   real(dp) function mult_binom_log_norm(n,m,s) result(ln)
      integer,intent(in)::n; real(dp),intent(in)::m,s
      real(dp)::mx,sumv,lw; integer::k
      mx=-huge(1.0_dp)
      do k=0,n; mx=max(mx,mult_binom_log_weight_cdf(k,n,m,s)); end do
      sumv=0.0_dp
      do k=0,n; lw=mult_binom_log_weight_cdf(k,n,m,s); sumv=sumv+exp(lw-mx); end do
      ln=mx+log(sumv)
   end function mult_binom_log_norm

   real(dp) function dmultbinom(y,size,m,s,log_p) result(v)
      ! The upstream C density and CDF use slightly different interaction
      ! expressions. This routine intentionally preserves the density path.
      integer,intent(in)::y,size; real(dp),intent(in)::m,s
      logical,intent(in),optional::log_p
      real(dp)::ss,lv; logical::lg
      lg=.false.; if(present(log_p))lg=log_p
      if(y<0 .or. y>size)then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
      ss=log(s)
      lv=log_choose(size,y)+real(y,dp)*log(m)+real(size-y,dp)* &
         (log(1.0_dp-m)+real((size-y)*y,dp)*ss)-mult_binom_log_norm(size,m,s)
      v=merge(lv,exp(lv),lg)
   end function dmultbinom

   real(dp) function pmultbinom(q,size,m,s) result(p)
      integer,intent(in)::q,size; real(dp),intent(in)::m,s
      real(dp)::ln; integer::k,qq
      if(q<0)then; p=0.0_dp; return; end if
      qq=min(q,size); ln=mult_binom_log_norm(size,m,s); p=0.0_dp
      do k=0,qq; p=p+exp(mult_binom_log_weight_cdf(k,size,m,s)-ln); end do
      p=min(1.0_dp,p)
   end function pmultbinom

   integer function qmultbinom(p,size,m,s) result(q)
      real(dp),intent(in)::p,m,s; integer,intent(in)::size
      q=discrete_quantile_bounded(p,size,cdf)
   contains
      function cdf(k) result(v); integer,intent(in)::k; real(dp)::v; v=pmultbinom(k,size,m,s); end function cdf
   end function qmultbinom

   function rmultbinom(n,size,m,s) result(x)
      integer,intent(in)::n,size; real(dp),intent(in)::m,s
      integer,allocatable::x(:); real(dp)::u; integer::i
      allocate(x(n)); do i=1,n; call random_number(u); x(i)=qmultbinom(u,size,m,s); end do
   end function rmultbinom

   real(dp) function double_pois_log_weight(y,m,s) result(lw)
      integer,intent(in)::y; real(dp),intent(in)::m,s
      integer::y1
      y1=max(y,1)
      lw=-s*m+real(y,dp)*s*(1.0_dp+log(m/real(y1,dp))) + &
         real(y,dp)*log(real(y1,dp))-real(y,dp)-log_gamma(real(y+1,dp))
   end function double_pois_log_weight

   real(dp) function double_pois_log_norm(my,m,s) result(ln)
      integer,intent(in)::my; real(dp),intent(in)::m,s
      real(dp)::mx,sumv,lw; integer::k
      mx=-huge(1.0_dp)
      do k=0,my; mx=max(mx,double_pois_log_weight(k,m,s)); end do
      sumv=0.0_dp
      do k=0,my; lw=double_pois_log_weight(k,m,s); sumv=sumv+exp(lw-mx); end do
      ln=mx+log(sumv)
   end function double_pois_log_norm

   real(dp) function ddoublepois(y,m,s,log_p) result(v)
      integer,intent(in)::y; real(dp),intent(in)::m,s
      logical,intent(in),optional::log_p
      integer::my; real(dp)::lv; logical::lg
      lg=.false.; if(present(log_p))lg=log_p
      if(y<0)then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
      my=3*max(y,100)
      lv=double_pois_log_weight(y,m,s)-double_pois_log_norm(my,m,s)
      v=merge(lv,exp(lv),lg)
   end function ddoublepois

   real(dp) function pdoublepois(q,m,s) result(p)
      integer,intent(in)::q; real(dp),intent(in)::m,s
      integer::my,k
      real(dp)::ln
      if(q<0)then; p=0.0_dp; return; end if
      my=3*max(q,100); ln=double_pois_log_norm(my,m,s); p=0.0_dp
      do k=0,q; p=p+exp(double_pois_log_weight(k,m,s)-ln); end do
      p=min(1.0_dp,p)
   end function pdoublepois

   integer function qdoublepois(p,m,s) result(q)
      real(dp),intent(in)::p,m,s
      q=discrete_quantile_unbounded(p,cdf)
   contains
      function cdf(k) result(v); integer,intent(in)::k; real(dp)::v; v=pdoublepois(k,m,s); end function cdf
   end function qdoublepois

   function rdoublepois(n,m,s) result(x)
      integer,intent(in)::n; real(dp),intent(in)::m,s
      integer,allocatable::x(:); real(dp)::u; integer::i
      allocate(x(n)); do i=1,n; call random_number(u); x(i)=qdoublepois(u,m,s); end do
   end function rdoublepois

   real(dp) function mult_pois_log_weight(y,m,s) result(lw)
      integer,intent(in)::y; real(dp),intent(in)::m,s
      lw=-m+real(y*y,dp)*log(s)+real(y,dp)*log(m)-log_gamma(real(y+1,dp))
   end function mult_pois_log_weight

   real(dp) function mult_pois_log_norm(my,m,s) result(ln)
      integer,intent(in)::my; real(dp),intent(in)::m,s
      real(dp)::mx,sumv,lw; integer::k
      mx=-huge(1.0_dp)
      do k=0,my; mx=max(mx,mult_pois_log_weight(k,m,s)); end do
      sumv=0.0_dp
      do k=0,my; lw=mult_pois_log_weight(k,m,s); sumv=sumv+exp(lw-mx); end do
      ln=mx+log(sumv)
   end function mult_pois_log_norm

   real(dp) function dmultpois(y,m,s,log_p) result(v)
      integer,intent(in)::y; real(dp),intent(in)::m,s
      logical,intent(in),optional::log_p
      integer::my; real(dp)::lv; logical::lg
      lg=.false.; if(present(log_p))lg=log_p
      if(y<0)then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
      my=3*max(y,100); lv=mult_pois_log_weight(y,m,s)-mult_pois_log_norm(my,m,s)
      v=merge(lv,exp(lv),lg)
   end function dmultpois

   real(dp) function pmultpois(q,m,s) result(p)
      integer,intent(in)::q; real(dp),intent(in)::m,s
      integer::my,k; real(dp)::ln
      if(q<0)then; p=0.0_dp; return; end if
      my=3*max(q,100); ln=mult_pois_log_norm(my,m,s); p=0.0_dp
      do k=0,q; p=p+exp(mult_pois_log_weight(k,m,s)-ln); end do
      p=min(1.0_dp,p)
   end function pmultpois

   integer function qmultpois(p,m,s) result(q)
      real(dp),intent(in)::p,m,s
      q=discrete_quantile_unbounded(p,cdf)
   contains
      function cdf(k) result(v); integer,intent(in)::k; real(dp)::v; v=pmultpois(k,m,s); end function cdf
   end function qmultpois

   function rmultpois(n,m,s) result(x)
      integer,intent(in)::n; real(dp),intent(in)::m,s
      integer,allocatable::x(:); real(dp)::u; integer::i
      allocate(x(n)); do i=1,n; call random_number(u); x(i)=qmultpois(u,m,s); end do
   end function rmultpois

   real(dp) function pvf_coeff_sum(y,m,s,f) result(r)
      integer,intent(in)::y; real(dp),intent(in)::m,s,f
      real(dp),allocatable::c(:,:)
      real(dp)::tmp1,tmp2,tmp3,tmp4
      integer::i,j
      if(y<=0)then; r=1.0_dp; return; end if
      allocate(c(0:y-1,0:y-1)); c=0.0_dp
      tmp1=gamma(1.0_dp-f); tmp2=log(m); tmp3=log(s+1.0_dp); tmp4=log(s)
      do i=0,y-1
         c(i,i)=1.0_dp
         if(i>0)then
            c(i,0)=gamma(real(i+1,dp)-f)/tmp1
            if(i>1)then
               do j=1,i-1
                  c(i,j)=c(i-1,j-1)+c(i-1,j)*(real(i,dp)-real(j+1,dp)*f)
               end do
            end if
         end if
      end do
      r=0.0_dp
      do i=1,y
         r=r+c(y-1,i-1)*exp(real(i,dp)*tmp2+(real(i,dp)*f-real(y,dp))*tmp3 - &
            real(i,dp)*(f-1.0_dp)*tmp4)
      end do
   end function pvf_coeff_sum

   real(dp) function dpvfpois(y,m,s,f,log_p) result(v)
      integer,intent(in)::y; real(dp),intent(in)::m,s,f
      logical,intent(in),optional::log_p
      real(dp)::pm,lv,size,prob; logical::lg
      lg=.false.; if(present(log_p))lg=log_p
      if(y<0)then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
      if(abs(f)<=epsilon(1.0_dp))then
         size=m*s; prob=s/(1.0_dp+s); pm=negative_binomial_pmf(y,size,prob)
      else
         pm=exp(-m*((s+1.0_dp)*((s+1.0_dp)/s)**(f-1.0_dp)-s)/f)
         if(y>0)pm=pm*pvf_coeff_sum(y,m,s,f)
         if(y>1)pm=pm/gamma(real(y+1,dp))
      end if
      if(lg)then; lv=log(max(pm,tiny(1.0_dp))); v=lv; else; v=pm; end if
   end function dpvfpois

   real(dp) function ppvfpois(q,m,s,f) result(p)
      integer,intent(in)::q; real(dp),intent(in)::m,s,f
      integer::k; real(dp)::size,prob
      if(q<0)then; p=0.0_dp; return; end if
      if(abs(f)<=epsilon(1.0_dp))then
         size=m*s; prob=s/(1.0_dp+s); p=negative_binomial_cdf(q,size,prob)
      else
         p=0.0_dp
         do k=0,q
            p=p+dpvfpois(k,m,s,f)
         end do
         p=min(1.0_dp,p)
      end if
   end function ppvfpois

   integer function qpvfpois(p,m,s,f) result(q)
      real(dp),intent(in)::p,m,s,f
      q=discrete_quantile_unbounded(p,cdf)
   contains
      function cdf(k) result(v); integer,intent(in)::k; real(dp)::v; v=ppvfpois(k,m,s,f); end function cdf
   end function qpvfpois

   function rpvfpois(n,m,s,f) result(x)
      integer,intent(in)::n; real(dp),intent(in)::m,s,f
      integer,allocatable::x(:); real(dp)::u; integer::i
      allocate(x(n)); do i=1,n; call random_number(u); x(i)=qpvfpois(u,m,s,f); end do
   end function rpvfpois

   real(dp) function pgammacount(q,m,s) result(p)
      integer,intent(in)::q; real(dp),intent(in)::m,s
      if(q<0)then; p=0.0_dp; else; p=regularized_gamma_q(real(q+1,dp)*s,m*s); end if
   end function pgammacount

   real(dp) function dgammacount(y,m,s,log_p) result(v)
      integer,intent(in)::y; real(dp),intent(in)::m,s
      logical,intent(in),optional::log_p
      real(dp)::pm; logical::lg
      lg=.false.; if(present(log_p))lg=log_p
      if(y<0)then; pm=0.0_dp
      else if(y==0)then; pm=regularized_gamma_q(s,m*s)
      else; pm=regularized_gamma_p(real(y,dp)*s,m*s)- &
         regularized_gamma_p(real(y+1,dp)*s,m*s)
      end if
      if(lg)then; v=log(max(pm,tiny(1.0_dp))); else; v=pm; end if
   end function dgammacount

   integer function qgammacount(p,m,s) result(q)
      real(dp),intent(in)::p,m,s
      q=discrete_quantile_unbounded(p,cdf)
   contains
      function cdf(k) result(v); integer,intent(in)::k; real(dp)::v; v=pgammacount(k,m,s); end function cdf
   end function qgammacount

   function rgammacount(n,m,s) result(x)
      integer,intent(in)::n; real(dp),intent(in)::m,s
      integer,allocatable::x(:); real(dp)::u; integer::i
      allocate(x(n)); do i=1,n; call random_number(u); x(i)=qgammacount(u,m,s); end do
   end function rgammacount

   real(dp) function dconsul(y,m,s,log_p) result(v)
      integer,intent(in)::y; real(dp),intent(in)::m,s
      logical,intent(in),optional::log_p
      real(dp)::base,lv; logical::lg
      lg=.false.; if(present(log_p))lg=log_p
      if(y<0)then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
      base=m+real(y,dp)*(s-1.0_dp)
      if(base<=0.0_dp)then; v=merge(-huge(1.0_dp),0.0_dp,lg); return; end if
      lv=log(m)-base/s+real(y-1,dp)*log(base)-real(y,dp)*log(s)-log_gamma(real(y+1,dp))
      v=merge(lv,exp(lv),lg)
   end function dconsul

   real(dp) function pconsul(q,m,s) result(p)
      integer,intent(in)::q; real(dp),intent(in)::m,s
      integer::k
      if(q<0)then; p=0.0_dp; return; end if
      p=0.0_dp; do k=0,q; p=p+dconsul(k,m,s); end do; p=min(1.0_dp,p)
   end function pconsul

   integer function qconsul(p,m,s) result(q)
      real(dp),intent(in)::p,m,s
      q=discrete_quantile_unbounded(p,cdf)
   contains
      function cdf(k) result(v); integer,intent(in)::k; real(dp)::v; v=pconsul(k,m,s); end function cdf
   end function qconsul

   function rconsul(n,m,s) result(x)
      integer,intent(in)::n; real(dp),intent(in)::m,s
      integer,allocatable::x(:); real(dp)::u; integer::i
      allocate(x(n)); do i=1,n; call random_number(u); x(i)=qconsul(u,m,s); end do
   end function rconsul

   integer function discrete_quantile_bounded(p,maxk,f) result(q)
      real(dp),intent(in)::p; integer,intent(in)::maxk
      interface
         function f(k) result(v)
            import dp
            integer,intent(in)::k
            real(dp)::v
         end function f
      end interface
      integer::lo,hi,mid
      if(p<=0.0_dp)then; q=0; return; else if(p>=1.0_dp)then; q=maxk; return; end if
      lo=0; hi=maxk
      do while(lo<hi)
         mid=(lo+hi)/2
         if(f(mid)>=p)then; hi=mid; else; lo=mid+1; end if
      end do
      q=lo
   end function discrete_quantile_bounded

   integer function discrete_quantile_unbounded(p,f) result(q)
      real(dp),intent(in)::p
      interface
         function f(k) result(v)
            import dp
            integer,intent(in)::k
            real(dp)::v
         end function f
      end interface
      integer::lo,hi,mid
      if(p<=0.0_dp)then; q=0; return; end if
      lo=0; hi=16
      do while(f(hi)<p .and. hi<1000000); hi=2*hi; end do
      do while(lo<hi)
         mid=(lo+hi)/2
         if(f(mid)>=p)then; hi=mid; else; lo=mid+1; end if
      end do
      q=lo
   end function discrete_quantile_unbounded

end module rmutil_discrete
