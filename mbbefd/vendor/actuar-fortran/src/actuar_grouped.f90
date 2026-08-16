! SPDX-License-Identifier: GPL-2.0-or-later
module actuar_grouped
  use actuar_kinds, only : dp
  use actuar_special, only : nan_dp
  implicit none
  private
  public :: grouped_mean, grouped_variance, grouped_quantile, ogive_cdf
  public :: empirical_moment, empirical_limited_moment, apply_coverage
contains

  pure function grouped_mean(breaks,frequency) result(m)
    real(dp), intent(in) :: breaks(:),frequency(:)
    real(dp) :: m
    real(dp), allocatable :: midpoint(:)
    if(size(breaks)/=size(frequency)+1 .or. sum(frequency)<=0.0_dp) then
      m=nan_dp(); return
    end if
    midpoint=0.5_dp*(breaks(:size(frequency))+breaks(2:))
    m=sum(midpoint*frequency)/sum(frequency)
  end function grouped_mean

  pure function grouped_variance(breaks,frequency,unbiased) result(v)
    real(dp), intent(in) :: breaks(:),frequency(:)
    logical, intent(in), optional :: unbiased
    real(dp) :: v,m,den
    real(dp), allocatable :: midpoint(:)
    logical :: ub
    if(size(breaks)/=size(frequency)+1 .or. sum(frequency)<=0.0_dp) then
      v=nan_dp(); return
    end if
    midpoint=0.5_dp*(breaks(:size(frequency))+breaks(2:))
    m=sum(midpoint*frequency)/sum(frequency)
    ub=.false.; if(present(unbiased)) ub=unbiased
    den=sum(frequency)-merge(1.0_dp,0.0_dp,ub)
    if(den<=0.0_dp) then
      v=nan_dp()
    else
      v=sum(frequency*(midpoint-m)**2)/den
    end if
  end function grouped_variance

  pure function ogive_cdf(x,breaks,frequency) result(p)
    real(dp), intent(in) :: x,breaks(:),frequency(:)
    real(dp) :: p,total,cum,fraction
    integer :: i
    if(size(breaks)/=size(frequency)+1 .or. sum(frequency)<=0.0_dp) then
      p=nan_dp(); return
    end if
    if(x<=breaks(1)) then; p=0.0_dp; return; end if
    if(x>=breaks(size(breaks))) then; p=1.0_dp; return; end if
    total=sum(frequency); cum=0.0_dp
    do i=1,size(frequency)
      if(x>=breaks(i+1)) then
        cum=cum+frequency(i)
      else if(x>breaks(i)) then
        fraction=(x-breaks(i))/(breaks(i+1)-breaks(i))
        cum=cum+frequency(i)*fraction
        exit
      end if
    end do
    p=cum/total
  end function ogive_cdf

  pure function grouped_quantile(prob,breaks,frequency) result(x)
    real(dp), intent(in) :: prob,breaks(:),frequency(:)
    real(dp) :: x,target,cum,prev,total
    integer :: i
    if(size(breaks)/=size(frequency)+1 .or. sum(frequency)<=0.0_dp .or. &
       prob<0.0_dp .or. prob>1.0_dp) then
      x=nan_dp(); return
    end if
    if(prob==0.0_dp) then; x=breaks(1); return; end if
    if(prob==1.0_dp) then; x=breaks(size(breaks)); return; end if
    total=sum(frequency); target=prob*total; cum=0.0_dp
    do i=1,size(frequency)
      prev=cum; cum=cum+frequency(i)
      if(cum>=target) then
        if(frequency(i)<=0.0_dp) then
          x=breaks(i)
        else
          x=breaks(i)+(breaks(i+1)-breaks(i))*(target-prev)/frequency(i)
        end if
        return
      end if
    end do
    x=breaks(size(breaks))
  end function grouped_quantile

  pure function empirical_moment(x,order,weights) result(m)
    real(dp), intent(in) :: x(:),order
    real(dp), intent(in), optional :: weights(:)
    real(dp) :: m
    if(size(x)==0) then; m=nan_dp(); return; end if
    if(present(weights)) then
      if(size(weights)/=size(x) .or. sum(weights)<=0.0_dp) then
        m=nan_dp()
      else
        m=sum(weights*x**order)/sum(weights)
      end if
    else
      m=sum(x**order)/real(size(x),dp)
    end if
  end function empirical_moment

  pure function empirical_limited_moment(x,limit,order,weights) result(m)
    real(dp), intent(in) :: x(:),limit,order
    real(dp), intent(in), optional :: weights(:)
    real(dp) :: m
    real(dp), allocatable :: y(:)
    allocate(y(size(x))); y=min(max(x,0.0_dp),limit)**order
    if(present(weights)) then
      if(size(weights)/=size(x) .or. sum(weights)<=0.0_dp) then
        m=nan_dp()
      else
        m=sum(weights*y)/sum(weights)
      end if
    else
      m=sum(y)/real(size(x),dp)
    end if
  end function empirical_limited_moment

  function apply_coverage(loss,deductible,limit,coinsurance,inflation,franchise,per_payment) result(payment)
    real(dp), intent(in) :: loss(:)
    real(dp), intent(in), optional :: deductible,limit,coinsurance,inflation
    logical, intent(in), optional :: franchise,per_payment
    real(dp), allocatable :: payment(:)
    real(dp) :: d,u,c,inf,x
    logical :: fr,pp
    integer :: i,nkeep
    real(dp), allocatable :: temp(:)
    d=0.0_dp; if(present(deductible)) d=deductible
    u=huge(1.0_dp); if(present(limit)) u=limit
    c=1.0_dp; if(present(coinsurance)) c=coinsurance
    inf=0.0_dp; if(present(inflation)) inf=inflation
    fr=.false.; if(present(franchise)) fr=franchise
    pp=.false.; if(present(per_payment)) pp=per_payment
    allocate(temp(size(loss))); temp=0.0_dp
    do i=1,size(loss)
      x=loss(i)*(1.0_dp+inf)
      if(fr) then
        if(x>d) temp(i)=c*min(x,u)
      else
        temp(i)=c*min(max(x-d,0.0_dp),max(0.0_dp,u-d))
      end if
    end do
    if(pp) then
      nkeep=count(temp>0.0_dp); allocate(payment(nkeep))
      payment=pack(temp,temp>0.0_dp)
    else
      allocate(payment(size(temp))); payment=temp
    end if
  end function apply_coverage

end module actuar_grouped
