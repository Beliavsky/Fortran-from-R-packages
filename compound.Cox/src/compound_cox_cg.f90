! SPDX-License-Identifier: GPL-2.0-only
module compound_cox_cg
  use compound_cox_kinds, only : dp
  use compound_cox_types, only : cg_result, cg_test_result
  use compound_cox_math, only : median_real, shuffle_int
  implicit none
  private
  public :: cg_clayton, cg_frank, cg_gumbel, cg_test
contains
  subroutine sort_pairs(time,status,ts,ds)
    real(dp),intent(in)::time(:)
    integer,intent(in)::status(:)
    real(dp),allocatable,intent(out)::ts(:)
    integer,allocatable,intent(out)::ds(:)
    integer::n,i,j
    real(dp)::tv
    integer::dv
    n=size(time)
    allocate(ts(n),ds(n))
    ts=time
    ds=status
    do i=2,n
    tv=ts(i)
    dv=ds(i)
    j=i-1
    do while(j>=1)
    if(ts(j)<=tv)exit
    ts(j+1)=ts(j)
    ds(j+1)=ds(j)
    j=j-1
    end do
    ts(j+1)=tv
    ds(j+1)=dv
    end do
  end subroutine sort_pairs

  subroutine finish_cg(ts, tau, surv, res)
    real(dp),intent(in)::ts(:),tau,surv(:)
    type(cg_result),intent(out)::res
    integer::i,n
    n=size(ts)
    allocate(res%time(n),res%n_risk(n),res%surv(n))
    res%time=ts
    res%surv=surv
    res%tau=tau
    do i=1,n
    res%n_risk(i)=real(n-i+1,dp)
    end do
    res%median=huge(1.0_dp)
    do i=1,n
    if(surv(i)<=0.5_dp)then
    res%median=ts(i)
    exit
    end if
    end do
  end subroutine finish_cg

  subroutine cg_clayton(time,status,alpha,res)
    real(dp),intent(in)::time(:),alpha
    integer,intent(in)::status(:)
    type(cg_result),intent(out)::res
    real(dp),allocatable::ts(:),a(:),s(:)
    integer,allocatable::ds(:)
    real(dp)::aa,cum
    integer::n,i,r
    aa=max(alpha,1.0e-9_dp)
    call sort_pairs(time,status,ts,ds)
    n=size(time)
    allocate(a(n),s(n))
    cum=0
    do i=1,n
      r=n-i+1
      if(i==n)then
      a(i)=0
      else
      a(i)=(real(r-1,dp)/n)**(-aa)-(real(r,dp)/n)**(-aa)
      end if
      cum=cum+a(i)*real(ds(i),dp)
      s(i)=(1.0_dp+cum)**(-1.0_dp/aa)
    end do
    call finish_cg(ts,aa/(aa+2.0_dp),s,res)
  end subroutine cg_clayton

  pure real(dp) function frank_integrand(x) result(v)
    real(dp),intent(in)::x
    if(abs(x)<1.0e-5_dp)then
    v=1.0_dp-x/2.0_dp+x*x/12.0_dp-x**4/720.0_dp
    else
    v=x/expm1_safe(x)
    end if
  end function frank_integrand
  pure real(dp) function expm1_safe(x) result(v)
    real(dp),intent(in)::x
    if(abs(x)<1e-5_dp)then
    v=x+x*x/2+x**3/6+x**4/24+x**5/120
    else
    v=exp(x)-1
    end if
  end function expm1_safe

  real(dp) function simpson_frank(a) result(val)
    real(dp),intent(in)::a
    integer,parameter::n=4096
    real(dp)::h,x,s
    integer::i
    if(abs(a)<1e-12_dp)then
    val=0
    return
    end if
    h=a/real(n,dp)
    s=frank_integrand(0.0_dp)+frank_integrand(a)
    do i=1,n-1
    x=h*real(i,dp)
    if(mod(i,2)==0)then
    s=s+2*frank_integrand(x)
    else
    s=s+4*frank_integrand(x)
    end if
    end do
    val=s*h/3
  end function simpson_frank

  subroutine cg_frank(time,status,alpha,res)
    real(dp),intent(in)::time(:),alpha
    integer,intent(in)::status(:)
    type(cg_result),intent(out)::res
    real(dp),allocatable::ts(:),a(:),s(:)
    integer,allocatable::ds(:)
    real(dp)::aa,cum,tau,intv,num,den
    integer::n,i,r
    aa=alpha
    if(abs(aa)<1e-6_dp)aa=1e-6_dp
    call sort_pairs(time,status,ts,ds)
    n=size(time)
    allocate(a(n),s(n))
    cum=0
    do i=1,n
      r=n-i+1
      if(i==n)then
      a(i)=0
      else
        num=exp(-aa*real(r-1,dp)/n)-1
        den=exp(-aa*real(r,dp)/n)-1
        a(i)=log(num/den)
      end if
      cum=cum+a(i)*real(ds(i),dp)
      s(i)=-log(1.0_dp+(exp(-aa)-1.0_dp)*exp(cum))/aa
    end do
    intv=simpson_frank(aa)
    tau=1.0_dp-4.0_dp/aa*(1.0_dp-intv/aa)
    call finish_cg(ts,tau,s,res)
  end subroutine cg_frank

  subroutine cg_gumbel(time,status,alpha,res)
    real(dp),intent(in)::time(:),alpha
    integer,intent(in)::status(:)
    type(cg_result),intent(out)::res
    real(dp),allocatable::ts(:),a(:),s(:)
    integer,allocatable::ds(:)
    real(dp)::aa,cum
    integer::n,i,r
    aa=max(alpha,0.0_dp)
    call sort_pairs(time,status,ts,ds)
    n=size(time)
    allocate(a(n),s(n))
    cum=0
    do i=1,n
      r=n-i+1
      if(i==n)then
      a(i)=0
      else
      a(i)=(-log(real(r-1,dp)/n))**(aa+1)-(-log(real(r,dp)/n))**(aa+1)
      end if
      cum=cum+a(i)*real(ds(i),dp)
      s(i)=exp(-cum**(1.0_dp/(1.0_dp+aa)))
    end do
    call finish_cg(ts,aa/(aa+1.0_dp),s,res)
  end subroutine cg_gumbel

  subroutine cg_dispatch(time,status,alpha,copula,res)
    real(dp),intent(in)::time(:),alpha
    integer,intent(in)::status(:)
    character(len=*),intent(in)::copula
    type(cg_result),intent(out)::res
    select case(copula(1:1))
    case('c','C')
    call cg_clayton(time,status,alpha,res)
    case('f','F')
    call cg_frank(time,status,alpha,res)
    case default
    call cg_gumbel(time,status,alpha,res)
    end select
  end subroutine cg_dispatch

  subroutine cg_group_stats(time,status,pi,cutoff,alpha,copula,tau,rmst_good,rmst_poor,l1)
    real(dp),intent(in)::time(:),pi(:),cutoff,alpha
    integer,intent(in)::status(:)
    character(len=*),intent(in)::copula
    real(dp),intent(out)::tau,rmst_good,rmst_poor,l1
    integer::n,ng,np,i,ig,ip,m
    logical,allocatable::good(:)
    real(dp),allocatable::tg(:),tp(:),ts(:),ssg(:),ssp(:)
    integer,allocatable::dg(:),dpv(:),idx(:)
    type(cg_result)::rg,rp
    real(dp)::tau_end,dt
    n=size(time)
    allocate(good(n))
    good=pi<=cutoff
    ng=count(good)
    np=n-ng
    allocate(tg(ng),dg(ng),tp(np),dpv(np))
    ig=0
    ip=0
    do i=1,n
    if(good(i))then
    ig=ig+1
    tg(ig)=time(i)
    dg(ig)=status(i)
    else
    ip=ip+1
    tp(ip)=time(i)
    dpv(ip)=status(i)
    end if
    end do
    call cg_dispatch(tg,dg,alpha,copula,rg)
    call cg_dispatch(tp,dpv,alpha,copula,rp)
    tau=rg%tau
    tau_end=min(maxval(tg),maxval(tp))
    allocate(idx(n))
    idx=[(i,i=1,n)]
    call sort_index(time,idx)
    m=count(time<=tau_end)
    allocate(ts(m),ssg(m),ssp(m))
    ig=0
    ip=0
    do i=1,m
      ts(i)=time(idx(i))
      if(good(idx(i)))then
      ig=ig+1
      else
      ip=ip+1
      end if
      if(ig>0)then
      ssg(i)=rg%surv(ig)
      else
      ssg(i)=1
      end if
      if(ip>0)then
      ssp(i)=rp%surv(ip)
      else
      ssp(i)=1
      end if
    end do
    rmst_good=0
    rmst_poor=0
    l1=0
    do i=1,m
      if(i==1)then
      dt=ts(1)
      rmst_good=rmst_good+dt
      rmst_poor=rmst_poor+dt
      else
      dt=ts(i)-ts(i-1)
      rmst_good=rmst_good+ssg(i-1)*dt
      rmst_poor=rmst_poor+ssp(i-1)*dt
      l1=l1+abs(ssg(i-1)-ssp(i-1))*dt
      end if
    end do
  end subroutine cg_group_stats

  subroutine sort_index(x,idx)
    real(dp),intent(in)::x(:)
    integer,intent(inout)::idx(:)
    integer::i,j,v
    do i=2,size(idx)
    v=idx(i)
    j=i-1
    do while(j>=1)
    if(x(idx(j))<=x(v))exit
    idx(j+1)=idx(j)
    j=j-1
    end do
    idx(j+1)=v
    end do
  end subroutine sort_index

  subroutine cg_test(time,status,pi,alpha,res,cutoff,copula,nperm)
    real(dp),intent(in)::time(:),pi(:),alpha
    integer,intent(in)::status(:)
    type(cg_test_result),intent(out)::res
    real(dp),intent(in),optional::cutoff
    character(len=*),intent(in),optional::copula
    integer,intent(in),optional::nperm
    real(dp)::cut,tau,rg,rp,l1,md,ml1,rgp,rpp,l1p,mdp,ml1p
    character(len=16)::cop
    integer::nper,j,n,cm,cl
    integer,allocatable::perm(:)
    real(dp),allocatable::pip(:)
    cut=median_real(pi)
    if(present(cutoff))cut=cutoff
    cop='clayton'
    if(present(copula))cop=copula
    nper=10000
    if(present(nperm))nper=nperm
    call cg_group_stats(time,status,pi,cut,alpha,cop,tau,rg,rp,l1)
    md=(rg-rp)/min_group_max(time,pi,cut)
    ml1=l1/min_group_max(time,pi,cut)
    cm=0
    cl=0
    n=size(time)
    allocate(perm(n),pip(n))
    perm=[(j,j=1,n)]
    do j=1,nper
    call shuffle_int(perm)
    pip=pi(perm)
    call cg_group_stats(time,status,pip,cut,alpha,cop,tau,rgp,rpp,l1p)
    mdp=(rgp-rpp)/min_group_max(time,pip,cut)
    ml1p=l1p/min_group_max(time,pip,cut)
    if(abs(mdp)>abs(md))cm=cm+1
    if(ml1p>ml1)cl=cl+1
    end do
    res%survival_diff=md
    res%rmstd=rg-rp
    res%p_value=real(cm,dp)/nper
    res%l1_distance=ml1
    res%integrated_l1=l1
    res%l1_p_value=real(cl,dp)/nper
    res%tau=tau
    res%n_good=count(pi<=cut)
    res%n_poor=n-res%n_good
    res%events_good=sum(status,mask=pi<=cut)
    res%events_poor=sum(status,mask=pi>cut)
    res%rmst_good=rg
    res%rmst_poor=rp
    res%mean_pi_good=sum(pi,mask=pi<=cut)/real(res%n_good,dp)
    res%mean_pi_poor=sum(pi,mask=pi>cut)/real(res%n_poor,dp)
  end subroutine cg_test

  real(dp) function min_group_max(time,pi,cut) result(tau)
    real(dp),intent(in)::time(:),pi(:),cut
    tau=min(maxval(time,mask=pi<=cut),maxval(time,mask=pi>cut))
  end function min_group_max
end module compound_cox_cg
