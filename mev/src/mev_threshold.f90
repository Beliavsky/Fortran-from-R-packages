module mev_threshold
  use mev_kinds, only: dp
  use mev_math, only: sort_descending
  use mev_tailindex, only: shape_hill
  use expint_mod, only: expint
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: bab_result, shape_lthill, shape_lthill_path, shape_trimhill, thselect_bab, bab_fcst

  type :: bab_result
    integer :: k0 = 0
    integer :: k0_lth = 0
    real(dp) :: thresh0 = 0.0_dp
    real(dp) :: shape = 0.0_dp
    real(dp) :: shape_lth = 0.0_dp
    real(dp) :: min_variance = 0.0_dp
    integer :: convergence = 1
    real(dp), allocatable :: stat(:)
    real(dp), allocatable :: ci_lower(:)
    real(dp), allocatable :: ci_upper(:)
    real(dp) :: test_size = 0.0_dp
    real(dp) :: test_level = 0.0_dp
    integer :: nsim = 0
  end type bab_result

contains

  real(dp) function shape_lthill(xdat,k,k0,sorted) result(xi)
    real(dp), intent(in) :: xdat(:)
    integer, intent(in) :: k,k0
    logical, intent(in), optional :: sorted
    real(dp), allocatable :: xs(:), tmp(:)
    real(dp) :: th, denom
    integer :: j
    logical :: srt
    xi=ieee_value(0.0_dp,ieee_quiet_nan)
    if(k<1 .or. k>=size(xdat) .or. k0<1 .or. k0>k) return
    srt=.false.; if(present(sorted)) srt=sorted
    allocate(xs(size(xdat)))
    if(srt) then
      xs=xdat
    else
      allocate(tmp(size(xdat))); call sort_descending(xdat,tmp); xs=tmp
    end if
    th=xs(k+1); if(th<=0.0_dp .or. any(xs(1:k0)<=0.0_dp)) return
    denom=1.0_dp
    if(k0<k) then
      do j=k0+1,k
        denom=denom+1.0_dp/real(j,dp)
      end do
    end if
    xi=(sum(log(xs(1:k0)))/real(k0,dp)-log(th))/denom
  end function shape_lthill

  subroutine shape_lthill_path(xdat,k,shape,k0min,sorted)
    real(dp), intent(in) :: xdat(:)
    integer, intent(in) :: k
    real(dp), allocatable, intent(out) :: shape(:)
    integer, intent(in), optional :: k0min
    logical, intent(in), optional :: sorted
    real(dp), allocatable :: xs(:), tmp(:)
    real(dp) :: th, csum, denom
    integer :: j,kmin
    logical :: srt
    kmin=1; if(present(k0min)) kmin=max(1,k0min)
    if(k<1 .or. k>=size(xdat) .or. kmin>k) then
      allocate(shape(0)); return
    end if
    srt=.false.; if(present(sorted)) srt=sorted
    allocate(xs(size(xdat)))
    if(srt) then
      xs=xdat
    else
      allocate(tmp(size(xdat))); call sort_descending(xdat,tmp); xs=tmp
    end if
    allocate(shape(k-kmin+1)); shape=ieee_value(0.0_dp,ieee_quiet_nan)
    th=xs(k+1); if(th<=0.0_dp .or. any(xs(1:k)<=0.0_dp)) return
    csum=0.0_dp
    do j=1,k
      csum=csum+log(xs(j))
      if(j>=kmin) then
        denom=1.0_dp+harmonic_range(j+1,k)
        shape(j-kmin+1)=(csum/real(j,dp)-log(th))/denom
      end if
    end do
  end subroutine shape_lthill_path

  real(dp) function shape_trimhill(xdat,k,k0,sorted) result(xi)
    real(dp), intent(in) :: xdat(:)
    integer, intent(in) :: k,k0
    logical, intent(in), optional :: sorted
    real(dp), allocatable :: xs(:), tmp(:), logy(:)
    real(dp) :: th
    logical :: srt
    xi=ieee_value(0.0_dp,ieee_quiet_nan)
    if(k<1 .or. k>=size(xdat) .or. k0<0 .or. k0>=k) return
    srt=.false.; if(present(sorted)) srt=sorted
    allocate(xs(size(xdat)))
    if(srt) then
      xs=xdat
    else
      allocate(tmp(size(xdat))); call sort_descending(xdat,tmp); xs=tmp
    end if
    th=xs(k+1); if(th<=0.0_dp .or. any(xs(k0+1:k)<=0.0_dp)) return
    allocate(logy(k-k0)); logy=log(xs(k0+1:k))-log(th)
    xi=real(k0+1,dp)/real(k-k0,dp)*logy(1)
    if(size(logy)>1) xi=xi+sum(logy(2:))/real(k-k0,dp)
  end function shape_trimhill

  real(dp) function bab_fcst(rho) result(v)
    real(dp), intent(in) :: rho
    real(dp) :: e1,e1mp,e1m2p,p
    p=rho
    if(p>=0.0_dp .or. abs(p)<1.0e-12_dp .or. abs(1.0_dp-p)<1.0e-12_dp) then
      v=ieee_value(0.0_dp,ieee_quiet_nan); return
    end if
    e1=expint(1.0_dp)
    e1mp=expint(1.0_dp-p)
    e1m2p=expint(1.0_dp-2.0_dp*p)
    v=(1.0_dp-exp(1.0_dp-2.0_dp*p)*(1.0_dp-2.0_dp*p)*e1m2p &
       -exp(2.0_dp-2.0_dp*p)*e1mp*e1mp)/(p*p*(1.0_dp-p)**2) &
      +2.0_dp/(p*p*(1.0_dp-p))*(exp(2.0_dp-p)*e1mp*e1-1.0_dp &
       +exp(1.0_dp-p)*(1.0_dp-p)*e1mp) &
      +(1.0_dp-exp(1.0_dp)*e1-exp(2.0_dp)*e1*e1)/(p*p)
  end function bab_fcst

  subroutine thselect_bab(xdat,res,kmin,kmax,rho,test,nsim,level)
    real(dp), intent(in) :: xdat(:)
    type(bab_result), intent(out) :: res
    integer, intent(in), optional :: kmin,kmax
    real(dp), intent(in), optional :: rho
    logical, intent(in), optional :: test
    integer, intent(in), optional :: nsim
    real(dp), intent(in), optional :: level
    real(dp), allocatable :: xs(:), path(:)
    real(dp) :: r,fc,varv,bestvar,meanv,lev
    integer :: lo,hi,k,idx,bestk,k0,ns
    logical :: do_test
    res=bab_result()
    if(size(xdat)<12) return
    lo=max(10,int(floor(0.2_dp*real(size(xdat),dp)))); if(present(kmin)) lo=kmin
    hi=size(xdat)-1; if(present(kmax)) hi=min(kmax,size(xdat)-1)
    r=-1.0_dp; if(present(rho)) r=rho
    if(lo<10 .or. hi<lo .or. r>=0.0_dp) return
    allocate(xs(size(xdat))); call sort_descending(xdat,xs)
    if(xs(hi+1)<=0.0_dp) return
    bestvar=huge(1.0_dp); bestk=0; meanv=0.0_dp
    do k=lo,hi
      call shape_lthill_path(xs,k,path,1,.true.)
      if(size(path)<2) cycle
      meanv=sum(path)/real(size(path),dp)
      varv=sum((path-meanv)**2)/real(size(path)-1,dp)
      if(varv<bestvar) then
        bestvar=varv; bestk=k; res%shape_lth=meanv
      end if
    end do
    if(bestk==0) return
    fc=bab_fcst(r); if(.not.(fc>0.0_dp)) return
    k0=floor(real(bestk,dp)*(0.502727_dp/((1.0_dp-r)**2*fc))**(-1.0_dp/(1.0_dp-2.0_dp*r)))
    k0=max(1,min(k0,size(xs)-1))
    idx=k0
    res%k0=k0; res%k0_lth=bestk; res%thresh0=xs(idx)
    res%shape=shape_hill(xs,k0); res%min_variance=bestvar; res%convergence=0
    do_test=.false.;if(present(test)) do_test=test
    if(do_test) then
      ns=999;if(present(nsim))ns=nsim
      lev=0.95_dp;if(present(level))lev=level
      if(ns>=19 .and. lev>0.0_dp .and. lev<1.0_dp .and. bestk>=3) then
        call bab_mc_test(xs,bestk,ns,lev,res)
      end if
    end if
  end subroutine thselect_bab

  subroutine bab_mc_test(xs,k,nsim,level,res)
    real(dp),intent(in)::xs(:),level
    integer,intent(in)::k,nsim
    type(bab_result),intent(inout)::res
    real(dp),allocatable::sh(:),tstat(:),om(:),omr(:),rmat(:,:),pp(:),lpr(:),cs(:)
    real(dp),allocatable::lo(:),hi(:),col(:)
    real(dp)::u,cov,bestdiff,diff
    integer::i,j,ktail,besttail,inside
    call shape_lthill_path(xs,k,sh,1,.true.)
    allocate(tstat(k-1),om(k-1),omr(k-1),rmat(nsim,k-1),pp(k+1),lpr(k),cs(k))
    tstat=sh(2:k)/sh(1:k-1)
    do j=1,k-1
      om(j)=real(j,dp)*(1.0_dp+harmonic_range(j+1,k))
    end do
    do j=1,k-2
      omr(j)=om(j)/om(j+1)
    end do
    omr(k-1)=1.0_dp
    do i=1,nsim
      pp=0.0_dp
      do j=1,k+1
        call random_number(u);u=max(u,tiny(1.0_dp))
        if(j==1) then;pp(j)=-log(u);else;pp(j)=pp(j-1)-log(u);end if
      end do
      lpr=log(pp(1:k)/pp(k+1));cs(1)=lpr(1)
      do j=2,k;cs(j)=cs(j-1)+lpr(j);end do
      do j=1,k-1
        rmat(i,j)=omr(j)*(1.0_dp+lpr(j+1)/cs(j))
      end do
    end do
    allocate(lo(k-1),hi(k-1),col(nsim))
    bestdiff=huge(1.0_dp);besttail=1
    do ktail=1,max(1,(nsim-1)/2)
      do j=1,k-1
        col=rmat(:,j);call sort_ascending_local(col)
        lo(j)=col(ktail);hi(j)=col(nsim+1-ktail)
      end do
      inside=0
      do i=1,nsim
        if(all(rmat(i,:)>=lo .and. rmat(i,:)<=hi)) inside=inside+1
      end do
      cov=real(inside,dp)/real(nsim,dp);diff=abs(cov-level)
      if(diff<bestdiff) then;bestdiff=diff;besttail=ktail;end if
      if(cov<level .and. ktail>1) exit
    end do
    do j=1,k-1
      col=rmat(:,j);call sort_ascending_local(col)
      lo(j)=col(besttail);hi(j)=col(nsim+1-besttail)
    end do
    allocate(res%stat(k-1),res%ci_lower(k-1),res%ci_upper(k-1))
    res%ci_lower=lo;res%ci_upper=hi
    res%stat=(tstat-lo)/max(hi-lo,epsilon(1.0_dp))
    res%test_size=real(count(tstat<lo .or. tstat>hi),dp)/real(k-1,dp)
    res%test_level=level;res%nsim=nsim
  end subroutine bab_mc_test

  subroutine sort_ascending_local(x)
    real(dp),intent(inout)::x(:)
    integer::i,j
    real(dp)::key
    do i=2,size(x)
      key=x(i);j=i-1
      do while(j>=1)
        if(x(j)<=key)exit
        x(j+1)=x(j);j=j-1
      end do
      x(j+1)=key
    end do
  end subroutine sort_ascending_local

  pure real(dp) function harmonic_range(a,b) result(h)
    integer, intent(in) :: a,b
    integer :: j
    h=0.0_dp
    if(a>b) return
    do j=a,b
      h=h+1.0_dp/real(j,dp)
    end do
  end function harmonic_range

end module mev_threshold
