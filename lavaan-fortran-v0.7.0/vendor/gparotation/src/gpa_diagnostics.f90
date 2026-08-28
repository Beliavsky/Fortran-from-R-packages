module gpa_diagnostics
  use gpa_kinds, only: dp
  use gpa_linalg, only: eye, inverse_matrix, logabsdet
  implicit none
  private
  type, public :: fit_stats
    real(dp) :: dof=0.0_dp, chisq=0.0_dp, srmr=0.0_dp, rmsea=0.0_dp
    real(dp) :: rmsea_l=0.0_dp, rmsea_u=0.0_dp, alpha=0.10_dp
    integer :: info=0
  end type fit_stats
  type, public :: simplicity_stats
    real(dp) :: hoffman=0.0_dp, gini=0.0_dp, bentler=0.0_dp
  end type simplicity_stats
  public :: residual_matrix, calc_fitstats, calc_fsi, calc_auc, calc_simplicity, calc_hyperplane
contains
  subroutine residual_matrix(rmat,l,phi,res)
    real(dp),intent(in)::rmat(:,:),l(:,:),phi(:,:)
    real(dp),intent(out)::res(size(rmat,1),size(rmat,2))
    real(dp)::rhat(size(rmat,1),size(rmat,2))
    integer::i
    rhat=matmul(matmul(l,phi),transpose(l))
    do i=1,size(rhat,1)
    rhat(i,i)=1.0_dp
    end do
    res=rmat-rhat
    do i=1,size(res,1)
    res(i,i)=0.0_dp
    end do
  end subroutine residual_matrix

  subroutine calc_fitstats(rmat,l,phi,nobs,stats,alpha)
    real(dp),intent(in)::rmat(:,:),l(:,:),phi(:,:)
    integer,intent(in)::nobs
    type(fit_stats),intent(out)::stats
    real(dp),intent(in),optional::alpha
    real(dp)::rhat(size(rmat,1),size(rmat,2)),res(size(rmat,1),size(rmat,2))
    real(dp)::ri(size(rmat,1),size(rmat,2)),fmin,ldh,ldr,tail,maxncp
    integer::v,k,i,j,info,nlower
    v=size(l,1)
    k=size(l,2)
    if(present(alpha)) stats%alpha=alpha
    stats%dof=real((v-k)**2-(v+k),dp)/2.0_dp
    rhat=matmul(matmul(l,phi),transpose(l))
    do i=1,v
    rhat(i,i)=1.0_dp
    end do
    res=rmat-rhat
    nlower=v*(v-1)/2
    stats%srmr=0.0_dp
    do i=2,v
    do j=1,i-1
    stats%srmr=stats%srmr+res(i,j)**2
    end do
    end do
    stats%srmr=sqrt(stats%srmr/real(max(1,nlower),dp))
    ldh=logabsdet(rhat,info=info)
    if(info/=0) then
    stats%info=info
    return
    end if
    ldr=logabsdet(rmat,info=info)
    if(info/=0) then
    stats%info=info
    return
    end if
    call inverse_matrix(rhat,ri,info)
    if(info/=0) then
    stats%info=info
    return
    end if
    fmin=ldh+trace(matmul(rmat,ri))-ldr-real(v,dp)
    stats%chisq=fmin*real(nobs-1,dp)
    if(stats%dof>0.0_dp .and. nobs>1) then
      stats%rmsea=sqrt(max(stats%chisq/(stats%dof*real(nobs,dp))-1.0_dp/real(nobs-1,dp),0.0_dp))
      tail=stats%alpha/2.0_dp
      maxncp=max(real(nobs,dp),stats%chisq)+2.0_dp*real(nobs,dp)
      stats%rmsea_l=rmsea_bound(stats%chisq,stats%dof,real(nobs,dp),1.0_dp-tail,maxncp)
      stats%rmsea_u=rmsea_bound(stats%chisq,stats%dof,real(nobs,dp),tail,maxncp)
    end if
  end subroutine calc_fitstats

  subroutine calc_fsi(l,fsi,mean_fsi,min_fsi)
    real(dp),intent(in)::l(:,:)
    real(dp),intent(out)::fsi(size(l,2)),mean_fsi,min_fsi
    real(dp)::s2,s4
    integer::j,p
    p=size(l,1)
    do j=1,size(l,2)
      s2=sum(l(:,j)**2)
      s4=sum(l(:,j)**4)
      if(s2<=tiny(1.0_dp)) then
      fsi(j)=0.0_dp
      else
      fsi(j)=(real(p,dp)*s4-s2*s2)/(real(p-1,dp)*s2*s2)
      end if
    end do
    mean_fsi=sum(fsi)/real(size(fsi),dp)
    min_fsi=minval(fsi)
  end subroutine calc_fsi

  subroutine calc_auc(l,auc,auc_adj,mean_auc,min_auc)
    real(dp),intent(in)::l(:,:)
    real(dp),intent(out)::auc(size(l,2)),auc_adj(size(l,2)),mean_auc,min_auc
    real(dp)::v(size(l,1)),s,c
    integer::j,i,p
    p=size(l,1)
    do j=1,size(l,2)
      v=l(:,j)**2
      call sort_desc(v)
      s=sum(v)
      c=0.0_dp
      auc(j)=0.0_dp
      if(s>tiny(1.0_dp)) then
        do i=1,p
        c=c+v(i)/s
        auc(j)=auc(j)+c
        end do
        auc(j)=auc(j)/real(p,dp)
      end if
    end do
    auc_adj=auc-0.5_dp
    mean_auc=sum(auc)/real(size(auc),dp)
    min_auc=minval(auc)
  end subroutine calc_auc

  function calc_simplicity(l) result(s)
    real(dp),intent(in)::l(:,:)
    type(simplicity_stats)::s
    real(dp)::l2(size(l,1),size(l,2)),a(size(l,1)),asum,g(size(l,2))
    real(dp)::tl2,tl4,cs2,den
    integer::p,k,i,j
    p=size(l,1)
    k=size(l,2)
    l2=l*l
    s%hoffman=1.0_dp-sum(l2*(1.0_dp-l2))/(real(p*k,dp)*0.25_dp)
    do j=1,k
      a=abs(l(:,j))
      call sort_asc(a)
      asum=sum(a)
      if(asum<=tiny(1.0_dp)) then
      g(j)=0.0_dp
      else
        g(j)=1.0_dp
        do i=1,p
        g(j)=g(j)-2.0_dp*(real(p-i,dp)+0.5_dp)*a(i)/(real(p,dp)*asum)
        end do
      end if
    end do
    s%gini=sum(g)/real(k,dp)
    tl2=sum(l2)
    tl4=sum(l2*l2)
    cs2=sum(sum(l2,dim=1)**2)
    den=tl2*tl2/real(p,dp)-tl4/real(p,dp)
    if(abs(den)>epsilon(1.0_dp)) s%bentler=(cs2/real(p,dp)-tl4/real(p,dp))/den
  end function calc_simplicity

  subroutine calc_hyperplane(l,cutoff,hp,total,maxhp,pct)
    real(dp),intent(in)::l(:,:),cutoff
    integer,intent(out)::hp(size(l,2)),total,maxhp
    real(dp),intent(out)::pct
    integer::j
    do j=1,size(l,2)
    hp(j)=count(abs(l(:,j))<cutoff)
    end do
    total=sum(hp)
    maxhp=size(l,1)*(size(l,2)-1)
    if(maxhp>0) pct=100.0_dp*real(total,dp)/real(maxhp,dp)
  end subroutine calc_hyperplane

  pure function trace(a) result(v)
    real(dp),intent(in)::a(:,:)
    real(dp)::v
    integer::i
    v=0.0_dp
    do i=1,min(size(a,1),size(a,2))
    v=v+a(i,i)
    end do
  end function trace

  function rmsea_bound(chisq,df,nobs,target,maxncp) result(v)
    real(dp),intent(in)::chisq,df,nobs,target,maxncp
    real(dp)::v,lo,hi,mid,flo,fhi,fm
    integer::it
    lo=0.0_dp
    hi=maxncp
    flo=ncx2_cdf(chisq,df,lo)-target
    fhi=ncx2_cdf(chisq,df,hi)-target
    if(flo*fhi>0.0_dp) then
    v=0.0_dp
    return
    end if
    do it=1,80
      mid=0.5_dp*(lo+hi)
      fm=ncx2_cdf(chisq,df,mid)-target
      if(abs(fm)<1.0e-10_dp) exit
      if(flo*fm<=0.0_dp) then
      hi=mid
      fhi=fm
      else
      lo=mid
      flo=fm
      end if
    end do
    v=sqrt(max(mid,0.0_dp)/(max(nobs-1.0_dp,1.0_dp)*df))
  end function rmsea_bound

  function ncx2_cdf(x,df,ncp) result(p)
    real(dp),intent(in)::x,df,ncp
    real(dp)::p,w,lam,term
    integer::j
    if(x<=0.0_dp) then
    p=0.0_dp
    return
    end if
    lam=ncp/2.0_dp
    w=exp(-lam)
    p=w*gamma_p(df/2.0_dp,x/2.0_dp)
    term=w
    do j=1,10000
      term=term*lam/real(j,dp)
      p=p+term*gamma_p(df/2.0_dp+real(j,dp),x/2.0_dp)
      if(term<1.0e-14_dp .and. j>lam+20.0_dp) exit
    end do
    p=min(1.0_dp,max(0.0_dp,p))
  end function ncx2_cdf

  function gamma_p(a,x) result(p)
    real(dp),intent(in)::a,x
    real(dp)::p,ap,del,sumv,b,c,d,h,an
    integer::n
    if(x<=0.0_dp) then
    p=0.0_dp
    return
    end if
    if(x<a+1.0_dp) then
      ap=a
      del=1.0_dp/a
      sumv=del
      do n=1,10000
        ap=ap+1.0_dp
        del=del*x/ap
        sumv=sumv+del
        if(abs(del)<abs(sumv)*1.0e-14_dp) exit
      end do
      p=sumv*exp(-x+a*log(x)-log_gamma(a))
    else
      b=x+1.0_dp-a
      c=1.0_dp/tiny(1.0_dp)
      d=1.0_dp/b
      h=d
      do n=1,10000
        an=-real(n,dp)*(real(n,dp)-a)
        b=b+2.0_dp
        d=an*d+b
        if(abs(d)<tiny(1.0_dp)) d=tiny(1.0_dp)
        c=b+an/c
        if(abs(c)<tiny(1.0_dp)) c=tiny(1.0_dp)
        d=1.0_dp/d
        del=d*c
        h=h*del
        if(abs(del-1.0_dp)<1.0e-14_dp) exit
      end do
      p=1.0_dp-exp(-x+a*log(x)-log_gamma(a))*h
    end if
    p=min(1.0_dp,max(0.0_dp,p))
  end function gamma_p

  subroutine sort_asc(x)
    real(dp),intent(inout)::x(:)
    real(dp)::t
    integer::i,j
    do i=2,size(x)
      t=x(i)
      j=i-1
      do while(j>=1)
        if(x(j)<=t) exit
        x(j+1)=x(j)
        j=j-1
      end do
      x(j+1)=t
    end do
  end subroutine sort_asc
  subroutine sort_desc(x)
    real(dp),intent(inout)::x(:)
    real(dp)::t
    integer::i,j
    do i=2,size(x)
      t=x(i)
      j=i-1
      do while(j>=1)
        if(x(j)>=t) exit
        x(j+1)=x(j)
        j=j-1
      end do
      x(j+1)=t
    end do
  end subroutine sort_desc
end module gpa_diagnostics
