! SPDX-License-Identifier: MIT
module gkwdist_core
   use gkwdist_kinds, only : dp
   use gkwdist_math, only : log1mexp, safe_pow, log_beta, beta_cdf, beta_quantile, &
      nan_dp, posinf_dp, neginf_dp, finite_dp, expm1_stable
   use gkwdist_rng, only : beta_rng
   use gkwdist_ad, only : ad2, ad_const, ad_var, ad_log_gamma, ad_log1mexp, log, exp, &
      operator(+), operator(-), operator(*), operator(/)
   implicit none
   private

   integer, parameter, public :: fam_gkw=1, fam_bkw=2, fam_kkw=3, fam_ekw=4
   integer, parameter, public :: fam_mc=5, fam_kw=6, fam_beta=7
   public :: family_from_name, family_npar, family_name
   public :: dgkw_scalar, pgkw_scalar, qgkw_scalar, rgkw_scalar
   public :: family_full_parameters, family_nll, family_derivatives

contains

   pure integer function family_from_name(name) result(fam)
      character(len=*), intent(in) :: name
      character(len=len(name)) :: s
      integer :: i, k
      s=name
      do i=1,len(s)
         k=iachar(s(i:i))
         if (k>=iachar('A') .and. k<=iachar('Z')) s(i:i)=achar(k+32)
      end do
      select case(trim(adjustl(s)))
      case('gkw'); fam=fam_gkw
      case('bkw'); fam=fam_bkw
      case('kkw'); fam=fam_kkw
      case('ekw'); fam=fam_ekw
      case('mc'); fam=fam_mc
      case('kw'); fam=fam_kw
      case('beta'); fam=fam_beta
      case default; fam=0
      end select
   end function family_from_name

   pure integer function family_npar(fam) result(n)
      integer, intent(in) :: fam
      select case(fam)
      case(fam_gkw); n=5
      case(fam_bkw,fam_kkw); n=4
      case(fam_ekw,fam_mc); n=3
      case(fam_kw,fam_beta); n=2
      case default; n=0
      end select
   end function family_npar

   pure function family_name(fam) result(name)
      integer, intent(in) :: fam
      character(len=4) :: name
      select case(fam)
      case(fam_gkw); name='gkw '
      case(fam_bkw); name='bkw '
      case(fam_kkw); name='kkw '
      case(fam_ekw); name='ekw '
      case(fam_mc); name='mc  '
      case(fam_kw); name='kw  '
      case(fam_beta); name='beta'
      case default; name='    '
      end select
   end function family_name

   pure logical function valid_full(a,b,g,d,l)
      real(dp), intent(in) :: a,b,g,d,l
      valid_full = finite_dp(a) .and. finite_dp(b) .and. finite_dp(g) .and. finite_dp(d) .and. finite_dp(l) .and. &
         a>0.0_dp .and. b>0.0_dp .and. g>0.0_dp .and. d>=0.0_dp .and. l>0.0_dp
   end function valid_full

   pure elemental function dgkw_scalar(x,a,b,g,d,l,log_prob) result(ans)
      real(dp), intent(in) :: x,a,b,g,d,l
      logical, intent(in), optional :: log_prob
      real(dp) :: ans, logdens, logu, logv, logw, logxa
      logical :: lp
      lp=.false.; if (present(log_prob)) lp=log_prob
      if (.not. valid_full(a,b,g,d,l)) then
         ans=nan_dp(); return
      end if
      if (.not. finite_dp(x) .or. x<=0.0_dp .or. x>=1.0_dp) then
         ans=merge(neginf_dp(),0.0_dp,lp); return
      end if
      logxa=a*log(x)
      logu=log1mexp(logxa)
      logv=log1mexp(b*logu)
      logw=log1mexp(l*logv)
      if (.not. finite_dp(logu) .or. .not. finite_dp(logv) .or. .not. finite_dp(logw)) then
         ans=merge(neginf_dp(),0.0_dp,lp); return
      end if
      logdens=log(l)+log(a)+log(b)-log_beta(g,d+1.0_dp)+(a-1.0_dp)*log(x)+(b-1.0_dp)*logu+ &
         (g*l-1.0_dp)*logv+d*logw
      if (lp) then
         ans=logdens
      else
         ans=exp(logdens)
      end if
   end function dgkw_scalar

   pure elemental function pgkw_scalar(q,a,b,g,d,l,lower_tail,log_p) result(ans)
      real(dp), intent(in) :: q,a,b,g,d,l
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: ans, logu, logv, y, p
      logical :: lower,lp
      lower=.true.; lp=.false.
      if (present(lower_tail)) lower=lower_tail
      if (present(log_p)) lp=log_p
      if (.not. valid_full(a,b,g,d,l)) then
         ans=nan_dp(); return
      end if
      if (.not. finite_dp(q) .or. q<=0.0_dp) then
         p=merge(0.0_dp,1.0_dp,lower)
      else if (q>=1.0_dp) then
         p=merge(1.0_dp,0.0_dp,lower)
      else
         logu=log1mexp(a*log(q))
         logv=log1mexp(b*logu)
         y=exp(l*logv)
         p=beta_cdf(y,g,d+1.0_dp,lower)
      end if
      if (lp) then
         if (p<=0.0_dp) then
            ans=neginf_dp()
         else
            ans=log(p)
         end if
      else
         ans=p
      end if
   end function pgkw_scalar

   pure elemental function qgkw_scalar(p,a,b,g,d,l,lower_tail,log_p) result(ans)
      real(dp), intent(in) :: p,a,b,g,d,l
      logical, intent(in), optional :: lower_tail,log_p
      real(dp) :: ans, pp, y, v, u
      logical :: lower,lp
      lower=.true.; lp=.false.
      if (present(lower_tail)) lower=lower_tail
      if (present(log_p)) lp=log_p
      if (.not. valid_full(a,b,g,d,l)) then
         ans=nan_dp(); return
      end if
      if (lp) then
         if (p>0.0_dp) then
            ans=nan_dp(); return
         end if
         if (lower) then
            pp=exp(p)
         else
            pp=-expm1_stable(p)
         end if
      else
         pp=merge(p,1.0_dp-p,lower)
      end if
      if (.not. finite_dp(pp) .or. pp<0.0_dp .or. pp>1.0_dp) then
         ans=nan_dp(); return
      else if (pp<=0.0_dp) then
         ans=0.0_dp; return
      else if (pp>=1.0_dp) then
         ans=1.0_dp; return
      end if
      y=beta_quantile(pp,g,d+1.0_dp)
      v=safe_pow(y,1.0_dp/l)
      u=safe_pow(max(0.0_dp,1.0_dp-v),1.0_dp/b)
      ans=safe_pow(max(0.0_dp,1.0_dp-u),1.0_dp/a)
      ans=max(0.0_dp,min(1.0_dp,ans))
   end function qgkw_scalar

   function rgkw_scalar(a,b,g,d,l) result(ans)
      real(dp), intent(in) :: a,b,g,d,l
      real(dp) :: ans, y, v, u
      if (.not. valid_full(a,b,g,d,l)) then
         ans=nan_dp(); return
      end if
      y=beta_rng(g,d+1.0_dp)
      if (y<=0.0_dp) then
         ans=0.0_dp
      else if (y>=1.0_dp) then
         ans=1.0_dp
      else
         v=safe_pow(y,1.0_dp/l)
         u=safe_pow(max(0.0_dp,1.0_dp-v),1.0_dp/b)
         ans=safe_pow(max(0.0_dp,1.0_dp-u),1.0_dp/a)
         ans=max(0.0_dp,min(1.0_dp,ans))
      end if
   end function rgkw_scalar

   pure subroutine family_full_parameters(fam,par,a,b,g,d,l,ok)
      integer, intent(in) :: fam
      real(dp), intent(in) :: par(:)
      real(dp), intent(out) :: a,b,g,d,l
      logical, intent(out) :: ok
      ok=size(par)>=family_npar(fam)
      a=1.0_dp; b=1.0_dp; g=1.0_dp; d=0.0_dp; l=1.0_dp
      if (.not.ok) return
      select case(fam)
      case(fam_gkw); a=par(1); b=par(2); g=par(3); d=par(4); l=par(5)
      case(fam_bkw); a=par(1); b=par(2); g=par(3); d=par(4)
      case(fam_kkw); a=par(1); b=par(2); d=par(3); l=par(4)
      case(fam_ekw); a=par(1); b=par(2); l=par(3)
      case(fam_mc); g=par(1); d=par(2); l=par(3)
      case(fam_kw); a=par(1); b=par(2)
      case(fam_beta); g=par(1); d=par(2)
      case default; ok=.false.
      end select
      if (ok) ok=valid_full(a,b,g,d,l)
   end subroutine family_full_parameters

   function family_nll(fam,par,data) result(nll)
      integer, intent(in) :: fam
      real(dp), intent(in) :: par(:),data(:)
      real(dp) :: nll,a,b,g,d,l,ld
      logical :: ok
      integer :: i
      call family_full_parameters(fam,par,a,b,g,d,l,ok)
      if (.not.ok .or. size(data)==0) then
         nll=posinf_dp(); return
      end if
      nll=0.0_dp
      do i=1,size(data)
         if (.not.finite_dp(data(i)) .or. data(i)<=0.0_dp .or. data(i)>=1.0_dp) then
            nll=posinf_dp(); return
         end if
         ld=dgkw_scalar(data(i),a,b,g,d,l,.true.)
         if (.not.finite_dp(ld)) then
            nll=posinf_dp(); return
         end if
         nll=nll-ld
      end do
   end function family_nll

   pure subroutine family_ad_parameters(fam,par,a,b,g,d,l,ok)
      integer, intent(in) :: fam
      real(dp), intent(in) :: par(:)
      type(ad2), intent(out) :: a,b,g,d,l
      logical, intent(out) :: ok
      real(dp) :: ar,br,gr,dr,lr
      call family_full_parameters(fam,par,ar,br,gr,dr,lr,ok)
      a=ad_const(ar); b=ad_const(br); g=ad_const(gr); d=ad_const(dr); l=ad_const(lr)
      if (.not.ok) return
      select case(fam)
      case(fam_gkw)
         a=ad_var(ar,1); b=ad_var(br,2); g=ad_var(gr,3); d=ad_var(dr,4); l=ad_var(lr,5)
      case(fam_bkw)
         a=ad_var(ar,1); b=ad_var(br,2); g=ad_var(gr,3); d=ad_var(dr,4)
      case(fam_kkw)
         a=ad_var(ar,1); b=ad_var(br,2); d=ad_var(dr,3); l=ad_var(lr,4)
      case(fam_ekw)
         a=ad_var(ar,1); b=ad_var(br,2); l=ad_var(lr,3)
      case(fam_mc)
         g=ad_var(gr,1); d=ad_var(dr,2); l=ad_var(lr,3)
      case(fam_kw)
         a=ad_var(ar,1); b=ad_var(br,2)
      case(fam_beta)
         g=ad_var(gr,1); d=ad_var(dr,2)
      end select
   end subroutine family_ad_parameters

   pure function ad_log_beta(a,b) result(c)
      type(ad2), intent(in) :: a,b
      type(ad2) :: c
      c=ad_log_gamma(a)+ad_log_gamma(b)-ad_log_gamma(a+b)
   end function ad_log_beta

   subroutine family_derivatives(fam,par,data,nll,gradient,hessian,ok)
      integer, intent(in) :: fam
      real(dp), intent(in) :: par(:),data(:)
      real(dp), intent(out) :: nll
      real(dp), intent(out) :: gradient(:),hessian(:,:)
      logical, intent(out), optional :: ok
      type(ad2) :: a,b,g,d,l,ll,logu,logv,logw
      logical :: valid
      integer :: i,npar
      npar=family_npar(fam)
      valid=npar>0 .and. size(par)>=npar .and. size(gradient)>=npar .and. &
         size(hessian,1)>=npar .and. size(hessian,2)>=npar .and. size(data)>0
      if (valid) call family_ad_parameters(fam,par,a,b,g,d,l,valid)
      if (valid) then
         do i=1,size(data)
            if (.not.finite_dp(data(i)) .or. data(i)<=0.0_dp .or. data(i)>=1.0_dp) then
               valid=.false.; exit
            end if
         end do
      end if
      gradient=nan_dp(); hessian=nan_dp()
      if (.not.valid) then
         nll=posinf_dp(); if(present(ok)) ok=.false.; return
      end if
      ll=real(size(data),dp)*(log(l)+log(a)+log(b)-ad_log_beta(g,d+1.0_dp))
      do i=1,size(data)
         logu=ad_log1mexp(a*log(data(i)))
         logv=ad_log1mexp(b*logu)
         logw=ad_log1mexp(l*logv)
         ll=ll+(a-1.0_dp)*log(data(i))+(b-1.0_dp)*logu+(g*l-1.0_dp)*logv+d*logw
      end do
      ll=-ll
      nll=ll%v
      gradient(1:npar)=ll%g(1:npar)
      hessian(1:npar,1:npar)=ll%h(1:npar,1:npar)
      if(present(ok)) ok=.true.
   end subroutine family_derivatives

end module gkwdist_core
