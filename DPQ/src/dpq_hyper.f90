! Hypergeometric approximations and utilities from DPQ.
! SPDX-License-Identifier: GPL-2.0-or-later
module dpq_hyper
   use r_compat, only: dp, phyper, dhyper, pbinom, dbinom, pbeta, r_lgamma, r_lchoose
   use dpq_core, only: prob_output, log1p_dp
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   implicit none
   private
   public :: phyper_appr_as152, phyper_ibeta, phyper1_molenaar, phyper2_molenaar
   public :: phyper_peizer, hyper2binom_p, supp_hyper, phyper_bin_molenaar
   public :: phyper_bin, phyper_all_bin, dhyper_bin_molenaar
   public :: lfastchoose, f05lchoose, bern, lgamma_asymp, phyper_r, pdhyper, phyper_r2

contains

   pure elemental real(dp) function pnorm_std(x) result(v)
      real(dp), intent(in) :: x
      v=0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function pnorm_std

   pure elemental real(dp) function phyper_appr_as152(q,m,n,k) result(v)
      real(dp), intent(in) :: q,m,n,k
      real(dp) :: mm,mean,sig
      mm=m+n
      mean=k*m/mm
      sig=sqrt(mean*n/mm*(mm-k)/(mm-1.0_dp))
      v=pnorm_std((q+0.5_dp-mean)/sig)
   end function phyper_appr_as152

   pure elemental real(dp) function phyper_ibeta(q,m,n,k) result(v)
      real(dp), intent(in) :: q,m,n,k
      real(dp) :: ntot,p,np,xi,dc,cc,lam
      ntot=m+n
      p=m/ntot
      np=k*p
      xi=(k+m-1.0_dp-2.0_dp*np)/(ntot-2.0_dp)
      dc=(ntot-k)*(1.0_dp-p)+np-1.0_dp
      cc=k*(k-1.0_dp)*p*(m-1.0_dp)/((ntot-1.0_dp)*dc)
      lam=(ntot-2.0_dp)**2*np*(ntot-k)*(1.0_dp-p) &
         /((ntot-1.0_dp)*dc*(k+m-1.0_dp-2.0_dp*np))
      v=pbeta(1.0_dp-xi,lam-q+cc,q-cc+1.0_dp)
   end function phyper_ibeta

   pure elemental real(dp) function phyper1_molenaar(q,m,n,k) result(v)
      real(dp), intent(in) :: q,m,n,k
      real(dp) :: ntot,z
      ntot=m+n
      z=2.0_dp/sqrt(ntot-1.0_dp)*(sqrt((q+1.0_dp)*(ntot-m-k+q+1.0_dp)) &
         -sqrt((k-q)*(m-q)))
      v=pnorm_std(z)
   end function phyper1_molenaar

   pure elemental real(dp) function phyper2_molenaar(q,m,n,k) result(v)
      real(dp), intent(in) :: q,m,n,k
      real(dp) :: ntot,z
      ntot=m+n
      z=2.0_dp/sqrt(ntot)*(sqrt((q+0.75_dp)*(ntot-m-k+q+0.75_dp)) &
         -sqrt(max(0.0_dp,(k-q-0.25_dp)*(m-q-0.25_dp))))
      v=pnorm_std(z)
   end function phyper2_molenaar

   pure elemental real(dp) function phyper_peizer(q,m,n,k) result(v)
      real(dp), intent(in) :: q,m,n,k
      real(dp) :: ntot,mm,nn,r,s,nd,a,b,c,d,ad,bd,cd,dd,l,z,den
      ntot=m+n
      nn=m
      mm=n
      r=k
      s=ntot-k
      nd=ntot-1.0_dp/6.0_dp
      a=q+0.5_dp
      ad=q+2.0_dp/3.0_dp
      b=m-q-0.5_dp
      bd=m-q-1.0_dp/3.0_dp
      c=k-q-0.5_dp
      cd=k-q-1.0_dp/3.0_dp
      d=ntot-m-k+q+0.5_dp
      dd=ntot-m-k+q+2.0_dp/3.0_dp
      if(min(a,b,c,d)<=0.0_dp)then
         v=phyper(q,nint(m),nint(n),nint(k))
         return
      end if
      l=a*log((a*ntot)/(nn*r))+b*log((b*ntot)/(nn*s)) &
         +c*log((c*ntot)/(mm*r))+d*log((d*ntot)/(mm*s))
      den=(mm+1.0_dp/6.0_dp)*(nn+1.0_dp/6.0_dp) &
         *(r+1.0_dp/6.0_dp)*(s+1.0_dp/6.0_dp)*ntot
      z=sign(1.0_dp,ad*dd-bd*cd)*sqrt(max(0.0_dp,2.0_dp*l*mm*nn*r*s*nd/den))
      v=pnorm_std(z)
   end function phyper_peizer

   pure elemental real(dp) function hyper2binom_p(x,m,n,k) result(v)
      real(dp), intent(in) :: x,m,n,k
      real(dp) :: ntot,p,nn
      ntot=m+n
      p=m/ntot
      nn=ntot-(k-1.0_dp)/2.0_dp
      v=(m-x/2.0_dp)/nn-k*(x-k*p-0.5_dp)/(6.0_dp*nn*nn)
      v=max(0.0_dp,min(1.0_dp,v))
   end function hyper2binom_p

   pure subroutine supp_hyper(m,n,k,lo,hi)
      integer, intent(in) :: m,n,k
      integer, intent(out) :: lo,hi
      lo=max(0,k-n)
      hi=min(k,m)
   end subroutine supp_hyper

   pure elemental real(dp) function phyper_bin_molenaar(q,m,n,k,variant,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,m,n,k
      integer, intent(in), optional :: variant
      logical, intent(in), optional :: lower_tail,log_p
      integer :: iv
      logical :: lt,lp
      real(dp) :: p0
      iv=1
      lt=.true.
      lp=.false.
      if(present(variant))iv=variant
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      select case(iv)
      case(1)
         p0=pbinom(q,nint(k),hyper2binom_p(q,m,n,k))
         v=prob_output(p0,lt,lp)
      case(2)
         p0=pbinom(q,nint(m),hyper2binom_p(q,k,n-k+m,m))
         v=prob_output(p0,lt,lp)
      case(3)
         p0=pbinom(m-1.0_dp-q,nint(m),hyper2binom_p(m-1.0_dp-q,m+n-k,k,m))
         v=prob_output(p0,.not.lt,lp)
      case default
         p0=pbinom(k-1.0_dp-q,nint(k),hyper2binom_p(k-1.0_dp-q,n,m,k))
         v=prob_output(p0,.not.lt,lp)
      end select
   end function phyper_bin_molenaar

   pure elemental real(dp) function phyper_bin(q,m,n,k,variant,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,m,n,k
      integer, intent(in), optional :: variant
      logical, intent(in), optional :: lower_tail,log_p
      integer :: iv
      logical :: lt,lp
      real(dp) :: p0,nt
      iv=1
      lt=.true.
      lp=.false.
      nt=m+n
      if(present(variant))iv=variant
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      select case(iv)
      case(1)
         p0=pbinom(q,nint(k),m/nt)
         v=prob_output(p0,lt,lp)
      case(2)
         p0=pbinom(q,nint(m),k/nt)
         v=prob_output(p0,lt,lp)
      case(3)
         p0=pbinom(m-1.0_dp-q,nint(m),(nt-k)/nt)
         v=prob_output(p0,.not.lt,lp)
      case default
         p0=pbinom(k-1.0_dp-q,nint(k),n/nt)
         v=prob_output(p0,.not.lt,lp)
      end select
   end function phyper_bin

   pure function phyper_all_bin(q,m,n,k,improved,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,m,n,k
      logical, intent(in), optional :: improved,lower_tail,log_p
      real(dp) :: v(4)
      logical :: imp,lt,lp
      integer :: j
      imp=.false.
      lt=.true.
      lp=.false.
      if(present(improved))imp=improved
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      do j=1,4
         if(imp)then
            v(j)=phyper_bin_molenaar(q,m,n,k,j,lt,lp)
         else
            v(j)=phyper_bin(q,m,n,k,j,lt,lp)
         end if
      end do
   end function phyper_all_bin

   pure elemental real(dp) function dhyper_bin_molenaar(x,m,n,k,log_p) result(v)
      real(dp), intent(in) :: x,m,n,k
      logical, intent(in), optional :: log_p
      logical :: lp
      real(dp) :: d
      lp=.false.
      if(present(log_p))lp=log_p
      d=dbinom(x,nint(k),hyper2binom_p(x,m,n,k),lp)
      v=d
   end function dhyper_bin_molenaar

   pure elemental real(dp) function lfastchoose(n,k) result(v)
      real(dp), intent(in) :: n,k
      v=r_lgamma(n+1.0_dp)-r_lgamma(k+1.0_dp)-r_lgamma(n-k+1.0_dp)
   end function lfastchoose

   pure elemental real(dp) function f05lchoose(n,k) result(v)
      real(dp), intent(in) :: n,k
      v=lfastchoose(anint(n),anint(k))
   end function f05lchoose

   pure real(dp) function bern(n) result(v)
      integer, intent(in) :: n
      real(dp), allocatable :: a(:)
      integer :: m,j
      if(n<0)then
         v=ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      if(n==1)then
         v=0.5_dp
         return
      end if
      if(n>1 .and. mod(n,2)==1)then
         v=0.0_dp
         return
      end if
      allocate(a(0:n))
      do m=0,n
         a(m)=1.0_dp/real(m+1,dp)
         do j=m,1,-1
            a(j-1)=real(j,dp)*(a(j-1)-a(j))
         end do
      end do
      v=a(0)
      if(n==1)v=-v
      deallocate(a)
   end function bern

   pure real(dp) function lgamma_asymp(x,n) result(v)
      real(dp), intent(in) :: x
      integer, intent(in) :: n
      real(dp) :: s,ix2,bs
      integer :: j
      s=(x-0.5_dp)*log(x)-x+0.5_dp*log(2.0_dp*acos(-1.0_dp))
      if(n<=0)then
         v=s
         return
      end if
      ix2=1.0_dp/(x*x)
      bs=0.0_dp
      do j=n,1,-1
         bs=bern(2*j)/(real(2*j*(2*j-1),dp))+bs*ix2
      end do
      v=s+bs/x
   end function lgamma_asymp

   pure elemental real(dp) function phyper_r(q,m,n,k,lower_tail,log_p) result(v)
      real(dp), intent(in) :: q,m,n,k
      logical, intent(in), optional :: lower_tail,log_p
      logical :: lt,lp
      real(dp)::p0
      lt=.true.
      lp=.false.
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      p0=phyper(q,nint(m),nint(n),nint(k))
      v=prob_output(p0,lt,lp)
   end function phyper_r

   pure elemental real(dp) function pdhyper(q,m,n,k,log_p) result(v)
      ! Ratio phyper(q,m,n,k) / dhyper(q,m,n,k), evaluated by the
      ! backwards recurrence used by DPQ::pdhyper.
      real(dp), intent(in) :: q,m,n,k
      logical, intent(in), optional :: log_p
      logical :: lp
      real(dp) :: qq, sumv, term
      lp=.false.
      if(present(log_p)) lp=log_p
      qq=floor(q+1.0e-7_dp)
      sumv=0.0_dp
      term=1.0_dp
      do while(qq>0.0_dp .and. term>=epsilon(1.0_dp)*sumv)
         term=term*(qq*(n-k+qq)/(k+1.0_dp-qq)/(m+1.0_dp-qq))
         sumv=sumv+term
         qq=qq-1.0_dp
      end do
      if(lp) then
      v=log1p_dp(sumv)
      else
      v=1.0_dp+sumv
      end if
   end function pdhyper

   pure elemental real(dp) function phyper_r2(q,m,n,k,lower_tail,log_p) result(v)
      ! R's phyperR2(): tail swapping plus pdhyper recurrence.
      real(dp), intent(in) :: q,m,n,k
      logical, intent(in), optional :: lower_tail,log_p
      logical :: lt,lp
      real(dp) :: qq,mm,nn,kk,tmp,d,pd,p0
      lt=.true.
      lp=.false.
      if(present(lower_tail))lt=lower_tail
      if(present(log_p))lp=log_p
      qq=floor(q+1.0e-7_dp)
      mm=anint(m)
      nn=anint(n)
      kk=anint(k)
      if(mm<0.0_dp .or. nn<0.0_dp .or. kk<0.0_dp .or. kk>mm+nn)then
         v=ieee_value(0.0_dp,ieee_quiet_nan)
         return
      end if
      if(qq*(mm+nn)>kk*mm)then
         tmp=nn
         nn=mm
         mm=tmp
         qq=kk-qq-1.0_dp
         lt=.not.lt
      end if
      if(qq<0.0_dp .or. qq<kk-nn)then
      v=prob_output(0.0_dp,lt,lp)
      return
      end if
      if(qq>=mm .or. qq>=kk)then
      v=prob_output(1.0_dp,lt,lp)
      return
      end if
      d=dhyper(qq,nint(mm),nint(nn),nint(kk))
      if(d<=0.0_dp)then
      v=prob_output(0.0_dp,lt,lp)
      return
      end if
      pd=pdhyper(qq,mm,nn,kk,.false.)
      p0=d*pd
      v=prob_output(min(1.0_dp,max(0.0_dp,p0)),lt,lp)
   end function phyper_r2

end module dpq_hyper
