! SPDX-License-Identifier: GPL-2.0-or-later
module tolerance_nonparametric
  use tolerance_kinds, only : dp
  use tolerance_types, only : tolerance_interval
  use tolerance_math, only : beta_i, binom_cdf, binom_pmf, binom_quantile, normal_quantile
  implicit none
  private
  public :: npbeta_tol_int, nptol_int, nptol_hm_options, nptol_ym_options, np_order, distfree_sample_size
  public :: distfree_alpha, distfree_confidence, distfree_content
contains

  function npbeta_tol_int(x,beta,side,known_lower,known_upper) result(out)
    real(dp),intent(in)::x(:)
    real(dp),intent(in),optional::beta,known_lower,known_upper
    integer,intent(in),optional::side
    type(tolerance_interval)::out
    real(dp),allocatable::z(:)
    real(dp)::be,ui,li,lo,up
    integer::n,ne,ne2,sd
    be=0.95_dp;if(present(beta))be=beta;sd=1;if(present(side))sd=side
    allocate(z(size(x)));z=x;call sort_local(z);n=size(z);lo=z(1);up=z(n)
    if(present(known_lower))lo=known_lower;if(present(known_upper))up=known_upper
    ui=be*real(n+1,dp);li=(real(n,dp)-ui)/2.0_dp
    ne=min(ceiling(be*real(n+1,dp)),n);ne2=max(floor(real(n-ne,dp)/2.0_dp),1)
    if(sd==1)then
      if(ui<=real(n,dp))up=z(ne)
      if(real(n,dp)-ui+1.0_dp>=1.0_dp)lo=z(max(n-ne+1,1))
    else
      if(ui+li<=real(n,dp))up=z(min(ne+ne2,n))
      if(li>=1.0_dp)lo=z(ne2)
    end if
    out%alpha=1.0_dp-be;out%p=be;out%estimate=0.0_dp;out%lower=lo;out%upper=up
  end function npbeta_tol_int

  function nptol_int(x,alpha,p,side,method,known_lower,known_upper) result(out)
    real(dp),intent(in)::x(:)
    real(dp),intent(in),optional::alpha,p,known_lower,known_upper
    integer,intent(in),optional::side
    character(len=*),intent(in),optional::method
    type(tolerance_interval)::out
    real(dp),allocatable::z(:)
    real(dp)::a,pp,lo,up,val,best,bval,gamma,piw,u1,u2,eps
    integer::n,sd,r,s,t,br,bs,i,j,hm_ind,diff,v1,v2,nmin
    character(len=8)::meth
    a=0.05_dp;if(present(alpha))a=alpha;pp=0.99_dp;if(present(p))pp=p
    sd=1;if(present(side))sd=side;meth='WILKS';if(present(method))meth=adjustl(method)
    allocate(z(size(x)));z=x;call sort_local(z);n=size(z);lo=z(1);up=z(n)
    if(present(known_lower))lo=known_lower;if(present(known_upper))up=known_upper
    select case(trim(meth))
    case('WILKS')
      if(sd==1)then
        r=binom_quantile(a,n,1.0_dp-pp);s=n-r+1
        if(r>=1)lo=z(r);if(s<=n)up=z(s)
      else
        best=huge(1.0_dp);br=0
        do r=1,max(0,(n-1)/2)
          val=(1.0_dp-beta_i(pp,real(n-2*r+1,dp),real(2*r,dp)))-(1.0_dp-a)
          if(val>0.0_dp .and. val<best)then;best=val;br=r;end if
        end do
        if(br>0)then;lo=z(br);up=z(n-br+1);end if
      end if
    case('WALD')
      if(sd==1)then
        r=binom_quantile(a,n,1.0_dp-pp);s=n-r+1
        if(r>=1)lo=z(r);if(s<=n)up=z(s)
      else
        best=huge(1.0_dp);br=0;bs=0
        do t=2,n
          do s=1,t-1
            val=(1.0_dp-beta_i(pp,real(t-s,dp),real(n-t+s+1,dp)))-(1.0_dp-a)
            if(val>0.0_dp .and. val<best)then;best=val;br=s;bs=t;end if
          end do
        end do
        if(br>0)then;lo=z(br);up=z(bs);end if
      end if
    case('HM')
      best=huge(1.0_dp);hm_ind=n
      do i=0,n
        val=binom_cdf(i,n,pp)-(1.0_dp-a)
        if(val>0.0_dp .and. val<best)then;best=val;hm_ind=i;end if
      end do
      diff=n-hm_ind
      if(sd==2)then
        if(diff>1)then
          v1=floor(real(diff,dp)/2.0_dp);v2=v1
          if(2*v1/=diff)v2=v1+1
          lo=z(max(1,v1));up=z(min(n,n-v2+1))
        end if
      else
        br=0
        do i=0,n
          val=(1.0_dp-binom_cdf(i-1,n,1.0_dp-pp))-(1.0_dp-a)
          if(val>0.0_dp)br=i
        end do
        if(br>0)lo=z(br)
        bs=n+1
        do i=1,n+1
          val=binom_cdf(i-1,n,pp)-(1.0_dp-a)
          if(val>0.0_dp)then;bs=i;exit;end if
        end do
        if(bs<=n)up=z(bs)
      end if
    case('YM')
      ! Young-Mathew interpolation/extrapolation.  Preserve the package's
      ! order-statistic logic while returning its OS-based interval.
      nmin=distfree_sample_size(a,pp,sd)
      if(sd==2)then
        call hm_two_indices(n,pp,a,br,bs)
        if(br>0)then;lo=z(br);up=z(bs);end if
      else
        gamma=1.0_dp-a
        r=binom_quantile(a,n,1.0_dp-pp);s=n-r+1
        if(n>=nmin .and. s>=1 .and. s<n)then
          piw=(gamma-binom_cdf(n-s-1,n,pp))/max(binom_pmf(n-s,n,pp),tiny(1.0_dp))
          lo=piw*z(s+1)+(1.0_dp-piw)*z(s)
        else if(n<nmin)then
          piw=-(gamma-binom_cdf(n-1,n,pp))/max(binom_pmf(n-1,n,pp),tiny(1.0_dp))
          lo=piw*z(min(2,n))+(1.0_dp-piw)*z(1)
        else if(r>=1)then;lo=z(r)
        end if
        if(n>=nmin .and. r>1 .and. r<=n)then
          piw=(gamma-binom_cdf(r-2,n,pp))/max(binom_pmf(r-1,n,pp),tiny(1.0_dp))
          up=piw*z(r)+(1.0_dp-piw)*z(r-1)
        else if(n<nmin)then
          piw=-(gamma-binom_cdf(n-1,n,pp))/max(binom_pmf(n-1,n,pp),tiny(1.0_dp))
          up=piw*z(max(1,n-1))+(1.0_dp-piw)*z(n)
        else if(s<=n)then;up=z(s)
        end if
        ! FOS roots are exposed separately by callers if desired; computing
        ! the roots here also exercises the beta interpolation used upstream.
        u1=beta_interp_root(1.0_dp-pp,a,n,.true.)
        u2=beta_interp_root(pp,a,n,.false.)
        eps=(real(n+1,dp)*u1-floor(real(n+1,dp)*u1)) + &
             (real(n+1,dp)*u2-floor(real(n+1,dp)*u2))
        if(eps<0.0_dp)up=up
      end if
    case default
      ! keep sample range
    end select
    out%alpha=a;out%p=pp;out%estimate=0.0_dp;out%lower=lo;out%upper=up
  end function nptol_int


  subroutine nptol_hm_options(x,alpha,p,lower,upper,noptions)
    ! Hahn-Meeker two-sided intervals.  For an odd number of excluded
    ! observations the R implementation returns two equally valid intervals.
    real(dp),intent(in)::x(:)
    real(dp),intent(in),optional::alpha,p
    real(dp),intent(out)::lower(2),upper(2)
    integer,intent(out)::noptions
    real(dp),allocatable::z(:)
    real(dp)::a,pp,val,best
    integer::n,i,hm_ind,diff,v1(2),v2(2),j
    a=0.05_dp;if(present(alpha))a=alpha;pp=0.99_dp;if(present(p))pp=p
    n=size(x);allocate(z(n));z=x;call sort_local(z);best=huge(1.0_dp);hm_ind=n
    do i=0,n
      val=binom_cdf(i,n,pp)-(1.0_dp-a)
      if(val>0.0_dp .and. val<best)then;best=val;hm_ind=i;end if
    end do
    diff=n-hm_ind;lower=z(1);upper=z(n)
    if(diff<=1)then;noptions=1;return;end if
    if(mod(diff,2)==0)then
      noptions=1;v1(1)=diff/2;v2(1)=diff/2
    else
      noptions=2;v1=[diff/2,(diff+1)/2];v2=[v1(1)+1,v1(2)-1]
    end if
    do j=1,noptions
      lower(j)=z(max(1,v1(j)));upper(j)=z(min(n,n-v2(j)+1))
    end do
  end subroutine nptol_hm_options


  subroutine nptol_ym_options(x,alpha,p,side,os,fos,has_fos)
    ! Young-Mathew interpolation/extrapolation.  For a one-sided request,
    ! the upstream routine returns both OS-based (Beran-Hall) and FOS-based
    ! (Hutson) rows.  For a two-sided request it returns one interpolated row.
    real(dp),intent(in)::x(:)
    real(dp),intent(in),optional::alpha,p
    integer,intent(in),optional::side
    type(tolerance_interval),intent(out)::os,fos
    logical,intent(out)::has_fos
    real(dp),allocatable::z(:),idxvals(:)
    real(dp)::a,pp,gamma,piw,u1,u2,e1,e2,g1,g2,ol(2),ou(2),width(4)
    real(dp)::hmlo(2),hmup(2),candlo(4),candup(4),pl(2),pu(2)
    integer::n,sd,nmin,lidx,uidx,nopt,i,best,r(2),sidx(2)

    a=0.05_dp;if(present(alpha))a=alpha
    pp=0.99_dp;if(present(p))pp=p
    sd=1;if(present(side))sd=side
    n=size(x);allocate(z(n));z=x;call sort_local(z);gamma=1.0_dp-a
    os%alpha=a;os%p=pp;os%estimate=0.0_dp;os%lower=z(1);os%upper=z(n)
    fos=os;has_fos=.false.

    if(sd==1)then
      has_fos=.true.;nmin=distfree_sample_size(a,pp,1)
      if(n>=nmin)then
        lidx=binom_quantile(a,n,1.0_dp-pp);uidx=n-lidx+1
        lidx=max(1,min(n,lidx));uidx=max(1,min(n,uidx))
        piw=(gamma-binom_cdf(n-lidx-1,n,pp))/max(binom_pmf(n-lidx,n,pp),tiny(1.0_dp))
        if(lidx==n)then;os%lower=z(lidx)
        else;os%lower=piw*z(lidx+1)+(1.0_dp-piw)*z(lidx);end if
        piw=(gamma-binom_cdf(uidx-2,n,pp))/max(binom_pmf(uidx-1,n,pp),tiny(1.0_dp))
        if(uidx==1)then;os%upper=z(uidx)
        else;os%upper=piw*z(uidx)+(1.0_dp-piw)*z(uidx-1);end if
        u1=beta_interp_root(1.0_dp-pp,a,n,.true.)
        u2=beta_interp_root(pp,a,n,.false.)
        e1=real(n+1,dp)*u1-floor(real(n+1,dp)*u1)
        e2=real(n+1,dp)*u2-floor(real(n+1,dp)*u2)
        if(lidx==n)then;fos%lower=z(lidx)
        else;fos%lower=(1.0_dp-e1)*z(lidx)+e1*z(lidx+1);end if
        if(uidx==1)then;fos%upper=z(uidx)
        else;fos%upper=(1.0_dp-e2)*z(uidx-1)+e2*z(uidx);end if
      else
        piw=-(gamma-binom_cdf(n-1,n,pp))/max(binom_pmf(n-1,n,pp),tiny(1.0_dp))
        os%lower=piw*z(min(2,n))+(1.0_dp-piw)*z(1)
        os%upper=piw*z(max(1,n-1))+(1.0_dp-piw)*z(n)
        u1=beta_interp_root(1.0_dp-pp,a,n,.true.)
        u2=beta_interp_root(pp,a,n,.false.)
        e1=-(real(n+1,dp)*u1-floor(real(n+1,dp)*u1))
        e2=-(real(n+1,dp)*u2-floor(real(n+1,dp)*u2))
        fos%lower=(1.0_dp-e1)*z(1)+e1*z(min(2,n))
        fos%upper=(1.0_dp-e2)*z(n)+e2*z(max(1,n-1))
      end if
      return
    end if

    ! The two-sided Young-Mathew routine starts from the Hahn-Meeker
    ! admissible order-statistic pair(s), then linearly interpolates in the
    ! achieved binomial confidence probability.
    allocate(idxvals(n));do i=1,n;idxvals(i)=real(i,dp);end do
    call nptol_hm_options(idxvals,a,pp,hmlo,hmup,nopt)
    do i=1,nopt
      r(i)=max(1,min(n,nint(hmlo(i))));sidx(i)=max(1,min(n,nint(hmup(i))))
    end do
    g1=binom_cdf(sidx(1)-r(1)-1,n,pp)
    g2=binom_cdf(sidx(1)-r(1)-2,n,pp)
    if(nopt==2)then
      do i=1,2
        pl(i)=linpred(z(r(i)),z(min(n,r(i)+1)),g1,g2,gamma)
        pu(i)=linpred(z(sidx(i)),z(max(1,sidx(i)-1)),g1,g2,gamma)
      end do
      candlo=[pl(1),pl(2),z(r(1)),z(r(2))]
      candup=[z(sidx(1)),z(sidx(2)),pu(1),pu(2)]
      width=candup-candlo;best=maxloc(width,dim=1)
      os%lower=candlo(best);os%upper=candup(best)
    else
      ol(1)=linpred(z(r(1)),z(min(n,r(1)+1)),g1,g2,gamma)
      ou(1)=linpred(z(sidx(1)),z(max(1,sidx(1)-1)),g1,g2,gamma)
      candlo(1:2)=[ol(1),z(r(1))];candup(1:2)=[z(sidx(1)),ou(1)]
      if(g1>=gamma)then
        width(1:2)=candup(1:2)-candlo(1:2);best=minloc(width(1:2),dim=1)
        os%lower=candlo(best);os%upper=candup(best)
      else
        os%lower=ol(1);os%upper=ou(1)
      end if
    end if
    fos=os
  end subroutine nptol_ym_options

  integer function np_order(m,alpha,p) result(n)
    integer,intent(in)::m
    real(dp),intent(in),optional::alpha,p
    real(dp)::a,pp,val
    integer::lo,hi,mid
    a=0.05_dp;if(present(alpha))a=alpha;pp=0.99_dp;if(present(p))pp=p
    lo=m;hi=max(m+1,2*m)
    do
      val=beta_i(1.0_dp-pp,real(m,dp),real(hi-m+1,dp))-(1.0_dp-a)
      if(val>=0.0_dp .or. hi>=100000*m)exit
      hi=min(2*hi,100000*m)
    end do
    do while(lo<hi)
      mid=(lo+hi)/2
      val=beta_i(1.0_dp-pp,real(m,dp),real(mid-m+1,dp))-(1.0_dp-a)
      if(val>=0.0_dp)then;hi=mid;else;lo=mid+1;end if
    end do
    n=lo
  end function np_order

  integer function distfree_sample_size(alpha,p,side) result(n)
    real(dp),intent(in)::alpha,p
    integer,intent(in),optional::side
    integer::sd
    real(dp)::v
    sd=1;if(present(side))sd=side
    if(sd==1)then
      n=max(1,ceiling(log(alpha)/log(p)))
    else
      n=2
      do
        v=real(n,dp)*p**(n-1)-real(n-1,dp)*p**n
        if(v<=alpha)exit
        n=n+1;if(n>100000000)exit
      end do
    end if
  end function distfree_sample_size

  real(dp) function distfree_alpha(n,p,side) result(alpha)
    integer,intent(in)::n
    real(dp),intent(in)::p
    integer,intent(in),optional::side
    integer::sd
    sd=1;if(present(side))sd=side
    if(sd==1)then;alpha=p**n;else;alpha=real(n,dp)*p**(n-1)-real(n-1,dp)*p**n;end if
  end function distfree_alpha


  real(dp) function distfree_confidence(n,p,side) result(confidence)
    integer,intent(in)::n
    real(dp),intent(in)::p
    integer,intent(in),optional::side
    confidence=1.0_dp-distfree_alpha(n,p,side)
  end function distfree_confidence

  real(dp) function distfree_content(n,alpha,side) result(p)
    integer,intent(in)::n
    real(dp),intent(in)::alpha
    integer,intent(in),optional::side
    integer::sd,it
    real(dp)::lo,hi,mid,v
    sd=1;if(present(side))sd=side
    if(sd==1)then;p=alpha**(1.0_dp/real(n,dp));return;end if
    lo=0.0_dp;hi=1.0_dp
    do it=1,100
      mid=0.5_dp*(lo+hi);v=real(n,dp)*mid**(n-1)-real(n-1,dp)*mid**n
      if(v>alpha)then;hi=mid;else;lo=mid;end if
    end do
    p=0.5_dp*(lo+hi)
  end function distfree_content

  subroutine sort_local(a)
    real(dp),intent(inout)::a(:)
    real(dp)::v
    integer::i,j
    do i=2,size(a);v=a(i);j=i-1;do while(j>=1);if(a(j)<=v)exit;a(j+1)=a(j);j=j-1;end do;a(j+1)=v;end do
  end subroutine sort_local

  subroutine hm_two_indices(n,p,a,r,s)
    integer,intent(in)::n;real(dp),intent(in)::p,a;integer,intent(out)::r,s
    integer::i,ind,diff
    real(dp)::v,best
    best=huge(1.0_dp);ind=n
    do i=0,n
      v=binom_cdf(i,n,p)-(1.0_dp-a)
      if(v>0.0_dp .and. v<best)then;best=v;ind=i;end if
    end do
    diff=n-ind
    if(diff<=1)then;r=1;s=n;else;r=max(1,diff/2);s=min(n,n-(diff-r)+1);end if
  end subroutine hm_two_indices


  pure real(dp) function linpred(x1,x2,g1,g2,target) result(v)
    real(dp),intent(in)::x1,x2,g1,g2,target
    if(abs(g2-g1)<=tiny(1.0_dp))then
      v=0.5_dp*(x1+x2)
    else
      v=x1+(target-g1)*(x2-x1)/(g2-g1)
    end if
  end function linpred

  real(dp) function beta_interp_root(u,a,n,lower_case) result(root)
    real(dp),intent(in)::u,a;integer,intent(in)::n;logical,intent(in)::lower_case
    real(dp)::lo,hi,mid,v;integer::it
    lo=1.0e-6_dp;hi=1.0_dp-1.0e-6_dp
    do it=1,100
      mid=0.5_dp*(lo+hi)
      v=beta_i(u,real(n+1,dp)*mid,real(n+1,dp)*(1.0_dp-mid))
      if(lower_case)v=v-(1.0_dp-a);if(.not.lower_case)v=v-a
      if(v>0.0_dp)then;hi=mid;else;lo=mid;end if
    end do
    root=0.5_dp*(lo+hi)
  end function beta_interp_root
end module tolerance_nonparametric
