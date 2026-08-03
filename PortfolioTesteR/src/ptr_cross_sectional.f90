! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
module ptr_cross_sectional
  use ptr_kinds, only : dp
  use ptr_utils, only : is_finite, nan_dp, finite_mean, finite_sd, rank_vector
  use ptr_indicators, only : panel_returns_simple
  implicit none
  private
  public :: calc_relative_strength_rank, calc_spread_indicators
  public :: calc_correlation_dispersion, calc_market_breadth
  public :: rank_within_groups, group_breadth, group_relative_indicators

contains

  subroutine calc_relative_strength_rank(indicator, out, method)
    real(dp),intent(in)::indicator(:,:)
    real(dp),allocatable,intent(out)::out(:,:)
    character(len=*),intent(in),optional::method
    character(len=16)::m
    real(dp),allocatable::r(:)
    real(dp)::mu,sd
    integer::t
    m='percentile';if(present(method))m=trim(method)
    allocate(out(size(indicator,1),size(indicator,2)));out=nan_dp()
    do t=1,size(indicator,1)
      select case(m)
      case('rank')
        call rank_vector(indicator(t,:),r,descending=.false.,normalize=.false.);out(t,:)=r
      case('zscore','z-score')
        mu=finite_mean(indicator(t,:));sd=finite_sd(indicator(t,:))
        if(is_finite(sd).and.sd>0.0_dp)then
          where(is_finite(indicator(t,:)))out(t,:)=(indicator(t,:)-mu)/sd
        end if
      case default
        call rank_vector(indicator(t,:),r,descending=.false.,normalize=.true.);out(t,:)=r
      end select
    end do
  end subroutine calc_relative_strength_rank

  subroutine calc_spread_indicators(indicator,top_n,bottom_n,spread)
    real(dp),intent(in)::indicator(:,:)
    integer,intent(in)::top_n,bottom_n
    real(dp),allocatable,intent(out)::spread(:)
    real(dp),allocatable::r(:)
    integer::t,j,ntop,nbot,n
    real(dp)::a,b
    allocate(spread(size(indicator,1)));spread=nan_dp()
    do t=1,size(indicator,1)
      call rank_vector(indicator(t,:),r,descending=.false.,normalize=.false.)
      n=count(is_finite(r));if(n==0)cycle
      ntop=0;nbot=0;a=0.0_dp;b=0.0_dp
      do j=1,size(r)
        if(.not.is_finite(r(j)))cycle
        if(r(j)>real(n-top_n,dp))then;a=a+indicator(t,j);ntop=ntop+1;end if
        if(r(j)<=real(bottom_n,dp))then;b=b+indicator(t,j);nbot=nbot+1;end if
      end do
      if(ntop>0.and.nbot>0)spread(t)=a/real(ntop,dp)-b/real(nbot,dp)
    end do
  end subroutine calc_spread_indicators

  subroutine calc_correlation_dispersion(prices,lookback,average_correlation,dispersion)
    real(dp),intent(in)::prices(:,:)
    integer,intent(in)::lookback
    real(dp),allocatable,intent(out)::average_correlation(:),dispersion(:)
    real(dp),allocatable::r(:,:)
    real(dp)::mx,my,sxy,sxx,syy,c,sumc,sumr
    integer::t,i,j,k,npair,nret
    call panel_returns_simple(prices,r)
    allocate(average_correlation(size(prices,1)),dispersion(size(prices,1)))
    average_correlation=nan_dp();dispersion=nan_dp()
    do t=lookback+1,size(prices,1)
      sumc=0.0_dp;npair=0
      do i=1,size(prices,2)-1
        do j=i+1,size(prices,2)
          mx=0.0_dp;my=0.0_dp;nret=0
          do k=t-lookback+1,t
            if(is_finite(r(k,i)).and.is_finite(r(k,j)))then
              mx=mx+r(k,i);my=my+r(k,j);nret=nret+1
            end if
          end do
          if(nret<2)cycle
          mx=mx/real(nret,dp);my=my/real(nret,dp);sxy=0.0_dp;sxx=0.0_dp;syy=0.0_dp
          do k=t-lookback+1,t
            if(is_finite(r(k,i)).and.is_finite(r(k,j)))then
              sxy=sxy+(r(k,i)-mx)*(r(k,j)-my)
              sxx=sxx+(r(k,i)-mx)**2;syy=syy+(r(k,j)-my)**2
            end if
          end do
          if(sxx>0.0_dp.and.syy>0.0_dp)then;c=sxy/sqrt(sxx*syy);sumc=sumc+c;npair=npair+1;end if
        end do
      end do
      if(npair>0)average_correlation(t)=sumc/real(npair,dp)
      sumr=0.0_dp;nret=0
      do j=1,size(prices,2)
        if(is_finite(r(t,j)))then;sumr=sumr+r(t,j);nret=nret+1;end if
      end do
      if(nret>1)then
        mx=sumr/real(nret,dp);sxx=0.0_dp
        do j=1,size(prices,2);if(is_finite(r(t,j)))sxx=sxx+(r(t,j)-mx)**2;end do
        dispersion(t)=sqrt(sxx/real(nret-1,dp))
      end if
    end do
  end subroutine calc_correlation_dispersion

  subroutine calc_market_breadth(condition,breadth,min_assets)
    real(dp),intent(in)::condition(:,:)
    real(dp),allocatable,intent(out)::breadth(:)
    integer,intent(in),optional::min_assets
    integer::t,n,minimum
    minimum=1;if(present(min_assets))minimum=min_assets
    allocate(breadth(size(condition,1)));breadth=nan_dp()
    do t=1,size(condition,1)
      n=count(is_finite(condition(t,:)))
      if(n>=minimum)breadth(t)=real(count(is_finite(condition(t,:)).and.condition(t,:)>0.0_dp),dp)/real(n,dp)
    end do
  end subroutine calc_market_breadth

  subroutine rank_within_groups(indicator,groups,out,normalize)
    real(dp),intent(in)::indicator(:,:)
    integer,intent(in)::groups(:)
    real(dp),allocatable,intent(out)::out(:,:)
    logical,intent(in),optional::normalize
    logical::norm
    real(dp),allocatable::v(:),r(:)
    integer::t,g,j,ng
    norm=.true.;if(present(normalize))norm=normalize
    allocate(out(size(indicator,1),size(indicator,2)));out=nan_dp()
    if(size(groups)/=size(indicator,2))return
    ng=maxval(groups);allocate(v(size(indicator,2)))
    do t=1,size(indicator,1)
      do g=1,ng
        v=nan_dp()
        do j=1,size(groups);if(groups(j)==g)v(j)=indicator(t,j);end do
        call rank_vector(v,r,descending=.false.,normalize=norm)
        do j=1,size(groups);if(groups(j)==g)out(t,j)=r(j);end do
      end do
    end do
  end subroutine rank_within_groups

  subroutine group_breadth(condition,groups,breadth)
    real(dp),intent(in)::condition(:,:)
    integer,intent(in)::groups(:)
    real(dp),allocatable,intent(out)::breadth(:,:)
    integer::t,g,j,n,yes,ng
    ng=maxval(groups);allocate(breadth(size(condition,1),ng));breadth=nan_dp()
    do t=1,size(condition,1)
      do g=1,ng
        n=0;yes=0
        do j=1,size(groups)
          if(groups(j)==g.and.is_finite(condition(t,j)))then
            n=n+1;if(condition(t,j)>0.0_dp)yes=yes+1
          end if
        end do
        if(n>0)breadth(t,g)=real(yes,dp)/real(n,dp)
      end do
    end do
  end subroutine group_breadth

  subroutine group_relative_indicators(indicator,groups,out)
    real(dp),intent(in)::indicator(:,:)
    integer,intent(in)::groups(:)
    real(dp),allocatable,intent(out)::out(:,:)
    real(dp)::mu
    integer::t,g,j,n,ng
    ng=maxval(groups);allocate(out(size(indicator,1),size(indicator,2)));out=nan_dp()
    do t=1,size(indicator,1)
      do g=1,ng
        mu=0.0_dp;n=0
        do j=1,size(groups)
          if(groups(j)==g.and.is_finite(indicator(t,j)))then;mu=mu+indicator(t,j);n=n+1;end if
        end do
        if(n==0)cycle
        mu=mu/real(n,dp)
        do j=1,size(groups)
          if(groups(j)==g.and.is_finite(indicator(t,j)))out(t,j)=indicator(t,j)-mu
        end do
      end do
    end do
  end subroutine group_relative_indicators

end module ptr_cross_sectional
