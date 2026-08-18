! SPDX-License-Identifier: GPL-2.0-or-later
module tolerance_support
  use tolerance_kinds, only : dp
  use tolerance_types, only : tolerance_interval
  use tolerance_normal, only : k_factor, normtol_int, bayesnormtol_int
  use tolerance_math, only : normal_cdf, normal_quantile, sample_mean, sample_sd
  implicit none
  private
  public :: k_table, norm_oc_content, norm_oc_alpha
  public :: norm_ss_fw, norm_ss_dir, norm_ss_ygzo, bonf_combine, bonf_tol_int

  abstract interface
    function one_sided_tol_proc(p,alpha) result(out)
      import :: dp, tolerance_interval
      real(dp),intent(in)::p,alpha
      type(tolerance_interval)::out
    end function one_sided_tol_proc
  end interface
contains

  subroutine k_table(n,alpha,p,side,method,m,k)
    integer,intent(in)::n(:)
    real(dp),intent(in)::alpha(:),p(:)
    integer,intent(in),optional::side,m
    character(len=*),intent(in),optional::method
    real(dp),intent(out)::k(size(n),size(alpha),size(p))
    integer::i,j,l,sd,mm
    character(len=8)::meth
    sd=1;if(present(side))sd=side;mm=50;if(present(m))mm=m;meth='HE';if(present(method))meth=adjustl(method)
    do i=1,size(n);do j=1,size(alpha);do l=1,size(p)
      k(i,j,l)=k_factor(n(i),alpha(j),p(l),sd,trim(meth),m=mm)
    end do;end do;end do
  end subroutine k_table

  real(dp) function norm_oc_content(k,n,alpha,side,method,m) result(p)
    real(dp),intent(in)::k,alpha
    integer,intent(in)::n
    integer,intent(in),optional::side,m
    character(len=*),intent(in),optional::method
    real(dp)::lo,hi,mid,v
    integer::it,sd,mm
    character(len=8)::meth
    sd=1;if(present(side))sd=side;mm=50;if(present(m))mm=m;meth='HE';if(present(method))meth=adjustl(method)
    lo=1.0e-10_dp;hi=1.0_dp-1.0e-10_dp
    do it=1,80
      mid=0.5_dp*(lo+hi);v=k_factor(n,alpha,mid,sd,trim(meth),m=mm)-k
      if(v>0.0_dp)then;hi=mid;else;lo=mid;end if
    end do
    p=0.5_dp*(lo+hi)
  end function norm_oc_content

  real(dp) function norm_oc_alpha(k,n,p,side,method,m) result(alpha)
    real(dp),intent(in)::k,p
    integer,intent(in)::n
    integer,intent(in),optional::side,m
    character(len=*),intent(in),optional::method
    real(dp)::lo,hi,mid,v
    integer::it,sd,mm
    character(len=8)::meth
    sd=1;if(present(side))sd=side;mm=50;if(present(m))mm=m;meth='HE';if(present(method))meth=adjustl(method)
    lo=1.0e-10_dp;hi=1.0_dp-1.0e-10_dp
    do it=1,80
      mid=0.5_dp*(lo+hi);v=k_factor(n,mid,p,sd,trim(meth),m=mm)-k
      ! K decreases as alpha grows.
      if(v>0.0_dp)then;lo=mid;else;hi=mid;end if
    end do
    alpha=0.5_dp*(lo+hi)
  end function norm_oc_alpha

  integer function norm_ss_fw(alpha,p,delta,pprime,side,m) result(nout)
    real(dp),intent(in)::alpha,p,delta,pprime
    integer,intent(in),optional::side,m
    integer::sd,mm,lo,hi,mid
    real(dp)::vlo,vhi,vm
    sd=1;if(present(side))sd=side;mm=50;if(present(m))mm=m
    lo=2;vlo=criterion(lo)
    hi=4;vhi=criterion(hi)
    do while(vlo*vhi>0.0_dp .and. hi<10000000)
      lo=hi;vlo=vhi;hi=2*hi;vhi=criterion(hi)
    end do
    if(vlo*vhi>0.0_dp)then;nout=hi;return;end if
    do while(hi-lo>1)
      mid=(lo+hi)/2;vm=criterion(mid)
      if(vlo*vm<=0.0_dp)then;hi=mid;vhi=vm;else;lo=mid;vlo=vm;end if
    end do
    if(abs(criterion(lo))<abs(criterion(hi)))then;nout=lo;else;nout=hi;end if
  contains
    real(dp) function criterion(n) result(v)
      integer,intent(in)::n
      if(sd==1)then
        v=k_factor(n,1.0_dp-delta,pprime,1,'HE',m=mm)-k_factor(n,alpha,p,1,'HE',m=mm)
      else
        v=k_factor(n,alpha,p,2,'HE',m=mm)-k_factor(n,1.0_dp-delta,pprime,2,'HE',m=mm)
      end if
    end function criterion
  end function norm_ss_fw

  integer function norm_ss_dir(mu,sigma,alpha,p,side,spec_lower,spec_upper,m) result(nout)
    real(dp),intent(in)::mu,sigma
    real(dp),intent(in),optional::alpha,p,spec_lower,spec_upper
    integer,intent(in),optional::side,m
    real(dp)::a,pp,kasy
    integer::sd,mm,n
    logical::ok
    a=0.05_dp;if(present(alpha))a=alpha;pp=0.99_dp;if(present(p))pp=p;sd=1;if(present(side))sd=side
    mm=50;if(present(m))mm=m
    if(sd==2)then
      if(.not.(present(spec_lower).and.present(spec_upper)))then;nout=huge(1);return;end if
      kasy=normal_quantile((1.0_dp+pp)/2.0_dp)
      if(mu-kasy*sigma<=spec_lower .or. mu+kasy*sigma>=spec_upper)then;nout=huge(1);return;end if
    end if
    do n=2,1000000
      if(sd==1)then
        if(present(spec_lower))then
          ok=mu-k_factor(n,a,pp,1,'OCT',m=mm)*sigma>=spec_lower
        else if(present(spec_upper))then
          ok=mu+k_factor(n,a,pp,1,'OCT',m=mm)*sigma<=spec_upper
        else;ok=.false.;end if
      else
        ok=mu-k_factor(n,a,pp,2,'OCT',m=mm)*sigma>=spec_lower .and. &
           mu+k_factor(n,a,pp,2,'OCT',m=mm)*sigma<=spec_upper
      end if
      if(ok)then;nout=n;return;end if
    end do
    nout=huge(1)
  end function norm_ss_dir

  integer function norm_ss_ygzo(x,alpha,p,side,spec_lower,spec_upper,delta,pprime,m, &
       mu0,sig20,m0,n0) result(nout)
    real(dp),intent(in)::x(:)
    real(dp),intent(in),optional::alpha,p,spec_lower,spec_upper,delta,pprime,mu0,sig20,m0,n0
    integer,intent(in),optional::side,m
    real(dp)::a,pp,de,ppr,mu,s,cont
    integer::sd,mm,n
    type(tolerance_interval)::ti
    a=0.05_dp;if(present(alpha))a=alpha;pp=0.99_dp;if(present(p))pp=p;sd=1;if(present(side))sd=side
    mm=50;if(present(m))mm=m;n=size(x);mu=sample_mean(x);s=sample_sd(x)
    if(present(mu0).and.present(sig20).and.present(m0).and.present(n0))then
      ti=bayesnormtol_int(mu,s,n,a,pp,sd,'HE',mu0,sig20,m0,n0);mu=ti%estimate;s=sqrt(sig20)
    else;ti=normtol_int(x,a,pp,sd,'EXACT');end if
    if(present(pprime))then;ppr=pprime
    else if(sd==2 .and. present(spec_lower).and.present(spec_upper))then
      ppr=normal_cdf((spec_upper-mu)/s)-normal_cdf((spec_lower-mu)/s)
      if(ppr<=pp .or. ppr>=1.0_dp)ppr=(1.0_dp+pp)/2.0_dp
    else if(present(spec_lower))then
      ppr=1.0_dp-normal_cdf((spec_lower-mu)/s);if(ppr<=pp .or. ppr>=1.0_dp)ppr=(1.0_dp+pp)/2.0_dp
    else if(present(spec_upper))then
      ppr=normal_cdf((spec_upper-mu)/s);if(ppr<=pp .or. ppr>=1.0_dp)ppr=(1.0_dp+pp)/2.0_dp
    else;ppr=(1.0_dp+pp)/2.0_dp;end if
    if(present(delta))then;de=delta
    else if(sd==1 .and. present(spec_lower))then
      cont=1.0_dp-normal_cdf((ti%lower-mu)/s);de=abs(cont-pp)/pp
    else if(sd==1)then
      cont=normal_cdf((ti%upper-mu)/s);de=abs(cont-pp)/pp
    else
      cont=normal_cdf((ti%upper-mu)/s)-normal_cdf((ti%lower-mu)/s);de=abs(cont-pp)/pp
    end if
    nout=norm_ss_fw(a,pp,max(min(de,1.0_dp-1.0e-8_dp),1.0e-8_dp),ppr,sd,mm)
  end function norm_ss_ygzo


  function bonf_tol_int(fn,p1,p2,alpha) result(out)
    procedure(one_sided_tol_proc)::fn
    real(dp),intent(in),optional::p1,p2,alpha
    type(tolerance_interval)::out,lower,upper
    real(dp)::q1,q2,a
    q1=0.005_dp;if(present(p1))q1=p1
    q2=0.005_dp;if(present(p2))q2=p2
    a=0.05_dp;if(present(alpha))a=alpha
    lower=fn(1.0_dp-q1,a);upper=fn(1.0_dp-q2,a)
    out=bonf_combine(lower,upper,q1,q2)
    out%alpha=a
  end function bonf_tol_int

  function bonf_combine(lower_one_sided,upper_one_sided,p1,p2) result(out)
    type(tolerance_interval),intent(in)::lower_one_sided,upper_one_sided
    real(dp),intent(in),optional::p1,p2
    type(tolerance_interval)::out
    real(dp)::q1,q2
    q1=0.005_dp;if(present(p1))q1=p1;q2=0.005_dp;if(present(p2))q2=p2
    out%alpha=lower_one_sided%alpha;out%p=1.0_dp-q1-q2;out%estimate=lower_one_sided%estimate
    out%lower=lower_one_sided%lower;out%upper=upper_one_sided%upper
  end function bonf_combine
end module tolerance_support
