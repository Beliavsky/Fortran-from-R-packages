! distr-fortran -- non-graphical QQ confidence-band calculations.
! Computational translation of qqbounds.R and internals-qqplot.R from distr.
! SPDX-License-Identifier: GPL-3.0-or-later
!
! This module depends on distr_ks (GPL-2.0-or-later).  It is therefore kept
! outside the LGPL umbrella module `distr`; import it explicitly when needed.
module distr_qq
   use distr_kinds, only : dp, nan_dp
   use distr_special, only : normal_quantile_std
   use distr_core, only : distribution_t, binomial_dist
   use distr_ks, only : p_ks2_asymptotic, p_kolmogorov2x
   implicit none
   private

   real(dp), parameter :: distr_resolution = 1.0e-6_dp
   real(dp), parameter :: root_tol = 1.0e-9_dp

   public :: qq_ks_critical, qq_pointwise_offsets, qqbounds

contains

   real(dp) function qq_ks_critical(alpha,n,exact) result(crit)
      ! Critical value used by upstream .q2kolmogorov.
      !
      ! NOTE: for exact=.false. this intentionally preserves upstream distr
      ! behavior, which solves 1 - pKS2(c) = alpha.  This differs from the
      ! usual alpha-CDF convention used by the exact branch.
      real(dp), intent(in) :: alpha
      integer, intent(in) :: n
      logical, intent(in), optional :: exact
      logical :: ex
      real(dp) :: lo,hi,mid,flo,fhi,fmid
      integer :: it

      if (n <= 0 .or. alpha <= 0.0_dp .or. alpha >= 1.0_dp) then
         crit = nan_dp()
         return
      end if
      ex = n < 100
      if (present(exact)) ex=exact

      if (ex) then
         lo=max(1.0e-12_dp,0.01_dp)
         hi=min(1.0_dp-1.0e-12_dp,3.0_dp*(1.0_dp-0.01_dp)/sqrt(real(n,dp)))
         flo=p_kolmogorov2x(lo,n)-alpha
         fhi=p_kolmogorov2x(hi,n)-alpha
         if (flo*fhi > 0.0_dp) then
            lo=1.0e-12_dp; hi=1.0_dp-1.0e-12_dp
            flo=p_kolmogorov2x(lo,n)-alpha
            fhi=p_kolmogorov2x(hi,n)-alpha
         end if
         if (flo*fhi > 0.0_dp) then
            crit=nan_dp(); return
         end if
         do it=1,100
            mid=0.5_dp*(lo+hi)
            fmid=p_kolmogorov2x(mid,n)-alpha
            if (abs(fmid) <= root_tol .or. hi-lo <= root_tol) exit
            if (flo*fmid <= 0.0_dp) then
               hi=mid; fhi=fmid
            else
               lo=mid; flo=fmid
            end if
         end do
         crit=mid*sqrt(real(n,dp))
      else
         lo=0.01_dp; hi=3.0_dp*(1.0_dp-0.01_dp)
         flo=1.0_dp-p_ks2_asymptotic(lo,1.0e-9_dp)-alpha
         fhi=1.0_dp-p_ks2_asymptotic(hi,1.0e-9_dp)-alpha
         if (flo*fhi > 0.0_dp) then
            lo=1.0e-8_dp; hi=8.0_dp
            flo=1.0_dp-p_ks2_asymptotic(lo,1.0e-9_dp)-alpha
            fhi=1.0_dp-p_ks2_asymptotic(hi,1.0e-9_dp)-alpha
         end if
         if (flo*fhi > 0.0_dp) then
            crit=nan_dp(); return
         end if
         do it=1,100
            mid=0.5_dp*(lo+hi)
            fmid=1.0_dp-p_ks2_asymptotic(mid,1.0e-9_dp)-alpha
            if (abs(fmid) <= root_tol .or. hi-lo <= root_tol) exit
            if (flo*fmid <= 0.0_dp) then
               hi=mid; fhi=fmid
            else
               lo=mid; flo=fmid
            end if
         end do
         crit=mid
      end if
   end function qq_ks_critical

   subroutine qq_pointwise_offsets(x,d,n,alpha,offsets,exact,nosym)
      ! Return offsets on the sqrt(n) scale, matching upstream .q2pw.
      real(dp), intent(in) :: x(:)
      type(distribution_t), intent(in) :: d
      integer, intent(in) :: n
      real(dp), intent(in) :: alpha
      real(dp), intent(out) :: offsets(size(x),2)
      logical, intent(in), optional :: exact,nosym
      logical :: ex,ns
      integer :: i
      real(dp) :: p,ld,ro,t,delta,nm

      if (n <= 0 .or. alpha <= 0.0_dp .or. alpha >= 1.0_dp) then
         offsets=nan_dp(); return
      end if
      ex=n<100; if (present(exact)) ex=exact
      ns=.false.; if (present(nosym)) ns=nosym

      do i=1,size(x)
         p=d%cdf(x(i))
         if (ex) then
            nm=sqrt(real(n,dp))*max(abs(x(i)),abs(real(n,dp)-x(i)))+1.0_dp
            if (ns) then
               call shortest_binom_offsets(x(i),p,d,n,alpha,nm,offsets(i,1),offsets(i,2))
            else
               t=solve_binom_width(x(i),p,d,n,alpha,0.0_dp,nm)
               if (t /= t) then
                  offsets(i,:)=nan_dp()
               else
                  offsets(i,1)=-t; offsets(i,2)=t
               end if
            end if
         else
            if (p<=0.0_dp .or. p>=1.0_dp) then
               offsets(i,:)=nan_dp()
               cycle
            end if
            ld=d%density(x(i),log_value=.true.)
            if (ld /= ld) then
               offsets(i,:)=nan_dp()
               cycle
            end if
            ro=exp(0.5_dp*(log(p)+log(1.0_dp-p))-ld)*normal_quantile_std(0.5_dp*(1.0_dp+alpha))
            offsets(i,1)=-ro; offsets(i,2)=ro
         end if
      end do
   end subroutine qq_pointwise_offsets

   subroutine qqbounds(x,d,alpha,n,bounds,with_pointwise,with_simultaneous, &
                       exact_pointwise,exact_simultaneous,nosym_pointwise)
      ! Non-graphical counterpart of exported R qqbounds().
      ! Columns: simultaneous left/right, pointwise left/right.
      real(dp), intent(in) :: x(:)
      type(distribution_t), intent(in) :: d
      real(dp), intent(in) :: alpha
      integer, intent(in) :: n
      real(dp), intent(out) :: bounds(size(x),4)
      logical, intent(in), optional :: with_pointwise,with_simultaneous
      logical, intent(in), optional :: exact_pointwise,exact_simultaneous,nosym_pointwise
      logical :: wp,ws,ep,es,ns
      real(dp) :: crit,pl,pr,sqrtn,ptail
      real(dp), allocatable :: off(:,:)
      integer :: i

      bounds=nan_dp()
      if (n<=0) return
      wp=.true.; ws=.true.; ep=n<100; es=n<100; ns=.false.
      if (present(with_pointwise)) wp=with_pointwise
      if (present(with_simultaneous)) ws=with_simultaneous
      if (present(exact_pointwise)) ep=exact_pointwise
      if (present(exact_simultaneous)) es=exact_simultaneous
      if (present(nosym_pointwise)) ns=nosym_pointwise
      sqrtn=sqrt(real(n,dp))

      if (ws) then
         crit=qq_ks_critical(alpha,n,es)
         if (crit==crit) then
            do i=1,size(x)
               pr=d%cdf(x(i))
               pl=d%cdf_left(x(i))
               bounds(i,1)=d%quantile(max(pl-crit/sqrtn,distr_resolution))
               ptail=max(1.0_dp-pr-crit/sqrtn,distr_resolution)
               bounds(i,2)=d%quantile(ptail,lower_tail=.false.)
            end do
         end if
      end if

      if (wp) then
         allocate(off(size(x),2))
         call qq_pointwise_offsets(x,d,n,alpha,off,ep,ns)
         bounds(:,3)=x+off(:,1)/sqrtn
         bounds(:,4)=x+off(:,2)/sqrtn
      end if
   end subroutine qqbounds

   real(dp) function binom_ci_value(t,pb,x,delta,d,n,alpha) result(v)
      real(dp), intent(in) :: t,pb,x,delta,alpha
      type(distribution_t), intent(in) :: d
      integer, intent(in) :: n
      type(distribution_t) :: bu,bl
      real(dp) :: pu,pl,npb,dr,pr,pleft
      integer :: k

      pu=min(1.0_dp,max(0.0_dp,d%cdf(x+(t+delta)/sqrt(real(n,dp)))))
      pl=min(1.0_dp,max(0.0_dp,d%cdf_left(x-(t-delta)/sqrt(real(n,dp)))))
      npb=real(n,dp)*pb
      k=int(floor(npb))
      dr=0.0_dp
      if (abs(npb-real(nint(npb),dp)) <= 16.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(npb))) then
         ! Upstream uses dbinom(npb, n, pmax(pu,1)); pmax makes the probability 1.
         if (nint(npb)==n) dr=1.0_dp
      end if
      bu=binomial_dist(n,pu)
      bl=binomial_dist(n,pl)
      pr=bu%cdf(real(k,dp),lower_tail=.false.)+dr
      pleft=bl%cdf(real(k,dp),lower_tail=.false.)
      v=pr-pleft-alpha
   end function binom_ci_value

   real(dp) function solve_binom_width(x,pb,d,n,alpha,delta,nm) result(root)
      real(dp), intent(in) :: x,pb,alpha,delta,nm
      type(distribution_t), intent(in) :: d
      integer, intent(in) :: n
      real(dp) :: lo,hi,flo,fhi,mid,fmid
      integer :: it,j

      root=nan_dp()
      lo=1.0e-12_dp
      flo=binom_ci_value(lo,pb,x,delta,d,n,alpha)
      hi=max(1.0_dp,min(nm,4.0_dp))
      fhi=binom_ci_value(hi,pb,x,delta,d,n,alpha)
      do j=1,60
         if (flo*fhi<=0.0_dp) exit
         if (hi>=nm) exit
         hi=min(nm,2.0_dp*hi)
         fhi=binom_ci_value(hi,pb,x,delta,d,n,alpha)
      end do
      if (flo*fhi>0.0_dp) return

      do it=1,100
         mid=0.5_dp*(lo+hi)
         fmid=binom_ci_value(mid,pb,x,delta,d,n,alpha)
         if (abs(fmid)<=root_tol .or. hi-lo<=root_tol) exit
         if (flo*fmid<=0.0_dp) then
            hi=mid; fhi=fmid
         else
            lo=mid; flo=fmid
         end if
      end do
      root=mid
   end function solve_binom_width

   subroutine shortest_binom_offsets(x,pb,d,n,alpha,nm,left,right)
      real(dp), intent(in) :: x,pb,alpha,nm
      type(distribution_t), intent(in) :: d
      integer, intent(in) :: n
      real(dp), intent(out) :: left,right
      real(dp), parameter :: gr=0.6180339887498948482_dp
      real(dp) :: a,b,c,e,fc,fe,tc,te,delta,t
      integer :: it

      a=-nm; b=nm
      c=b-gr*(b-a); e=a+gr*(b-a)
      tc=solve_binom_width(x,pb,d,n,alpha,c,nm)
      te=solve_binom_width(x,pb,d,n,alpha,e,nm)
      fc=merge(tc,huge(1.0_dp),tc==tc)
      fe=merge(te,huge(1.0_dp),te==te)
      do it=1,100
         if (abs(b-a)<=root_tol*max(1.0_dp,nm)) exit
         if (fc<fe) then
            b=e; e=c; fe=fc
            c=b-gr*(b-a)
            tc=solve_binom_width(x,pb,d,n,alpha,c,nm)
            fc=merge(tc,huge(1.0_dp),tc==tc)
         else
            a=c; c=e; fc=fe
            e=a+gr*(b-a)
            te=solve_binom_width(x,pb,d,n,alpha,e,nm)
            fe=merge(te,huge(1.0_dp),te==te)
         end if
      end do
      delta=0.5_dp*(a+b)
      t=solve_binom_width(x,pb,d,n,alpha,delta,nm)
      if (t/=t) then
         left=nan_dp(); right=nan_dp()
      else
         left=-t+delta; right=t+delta
      end if
   end subroutine shortest_binom_offsets

end module distr_qq
