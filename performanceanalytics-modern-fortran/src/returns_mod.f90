! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
module returns_mod
  use kinds_mod, only: dp
  use statistics_mod, only: mean_value, sd_value, median_value
  implicit none
  private
  public :: calculate_returns, cumulative_return, annualized_return, annualized_excess_return
  public :: excess_returns, relative_returns, centered_returns, geltner_unsmooth
  public :: loc_scale_robust, clean_boudt, portfolio_result, portfolio_returns
  public :: wealth_index, convert_return_frequency, level_from_returns
  type :: portfolio_result
    real(dp), allocatable :: returns(:)
    real(dp), allocatable :: wealth(:)
    real(dp), allocatable :: contributions(:,:)
    real(dp), allocatable :: begin_weights(:,:)
    real(dp), allocatable :: end_weights(:,:)
    real(dp), allocatable :: turnover(:)
  end type portfolio_result
contains
  subroutine calculate_returns(prices,method,r)
    real(dp),intent(in)::prices(:)
    character(len=*),intent(in)::method
    real(dp),allocatable,intent(out)::r(:)
    integer::n
    n=size(prices); allocate(r(max(0,n-1)))
    if(n<2) return
    select case(trim(adjustl(method)))
    case('log','compound','continuous')
      r=log(prices(2:n)/prices(1:n-1))
    case('difference','diff')
      r=prices(2:n)-prices(1:n-1)
    case default
      r=prices(2:n)/prices(1:n-1)-1.0_dp
    end select
  end subroutine calculate_returns

  pure real(dp) function cumulative_return(r,geometric) result(v)
    real(dp),intent(in)::r(:)
    logical,intent(in),optional::geometric
    logical::g
    g=.true.; if(present(geometric)) g=geometric
    if(g) then
      if(any(1.0_dp+r<=0.0_dp)) then; v=-1.0_dp; else; v=exp(sum(log(1.0_dp+r)))-1.0_dp; end if
    else
      v=sum(r)
    end if
  end function cumulative_return

  pure real(dp) function annualized_return(r,scale,geometric) result(v)
    real(dp),intent(in)::r(:),scale
    logical,intent(in),optional::geometric
    logical::g
    g=.true.; if(present(geometric)) g=geometric
    if(size(r)==0) then; v=0.0_dp; return; end if
    if(g) then
      if(any(1.0_dp+r<=0.0_dp)) then; v=-1.0_dp; else; v=exp(scale/real(size(r),dp)*sum(log(1.0_dp+r)))-1.0_dp; end if
    else
      v=mean_value(r)*scale
    end if
  end function annualized_return

  pure real(dp) function annualized_excess_return(r,rf,scale,geometric) result(v)
    real(dp),intent(in)::r(:),rf(:),scale
    logical,intent(in),optional::geometric
    real(dp),allocatable::e(:)
    integer::n
    n=min(size(r),size(rf)); allocate(e(n)); e=r(:n)-rf(:n)
    if(present(geometric)) then; v=annualized_return(e,scale,geometric); else; v=annualized_return(e,scale); end if
  end function annualized_excess_return

  pure subroutine excess_returns(r,rf,out)
    real(dp),intent(in)::r(:),rf(:)
    real(dp),intent(out)::out(:)
    integer::n
    n=min(size(r),size(rf)); out(:n)=r(:n)-rf(:n)
  end subroutine excess_returns

  pure subroutine relative_returns(r,benchmark,out,geometric)
    real(dp),intent(in)::r(:),benchmark(:)
    real(dp),intent(out)::out(:)
    logical,intent(in),optional::geometric
    logical::g
    integer::n
    n=min(size(r),size(benchmark)); g=.true.; if(present(geometric)) g=geometric
    if(g) then
      out(:n)=(1.0_dp+r(:n))/(1.0_dp+benchmark(:n))-1.0_dp
    else
      out(:n)=r(:n)-benchmark(:n)
    end if
  end subroutine relative_returns

  pure subroutine centered_returns(r,out)
    real(dp),intent(in)::r(:)
    real(dp),intent(out)::out(:)
    out=r-mean_value(r)
  end subroutine centered_returns

  pure subroutine wealth_index(r,value,out)
    real(dp),intent(in)::r(:),value
    real(dp),intent(out)::out(:)
    integer::i
    if(size(out)==0) return
    out(1)=value*(1.0_dp+r(1))
    do i=2,min(size(r),size(out)); out(i)=out(i-1)*(1.0_dp+r(i)); end do
  end subroutine wealth_index

  pure subroutine geltner_unsmooth(r,out)
    real(dp),intent(in)::r(:)
    real(dp),intent(out)::out(:)
    real(dp)::rho
    integer::n
    n=size(r); if(n==0) return
    if(n<3) then; out=r; return; end if
    rho=sum((r(2:n)-mean_value(r(2:n)))*(r(1:n-1)-mean_value(r(1:n-1))))/ &
      max(sum((r(1:n-1)-mean_value(r(1:n-1)))**2),tiny(1.0_dp))
    rho=max(-0.99_dp,min(0.99_dp,rho)); out(1)=r(1)
    out(2:n)=(r(2:n)-rho*r(1:n-1))/(1.0_dp-rho)
  end subroutine geltner_unsmooth

  subroutine loc_scale_robust(r,location,scale)
    real(dp),intent(in)::r(:)
    real(dp),intent(out)::location,scale
    real(dp),allocatable::d(:)
    location=median_value(r); allocate(d(size(r))); d=abs(r-location)
    scale=1.4826_dp*median_value(d)
    if(scale<=tiny(1.0_dp)) scale=sd_value(r)
  end subroutine loc_scale_robust

  subroutine clean_boudt(r,alpha,out)
    real(dp),intent(in)::r(:),alpha
    real(dp),intent(out)::out(:)
    real(dp)::location,scale,cut
    call loc_scale_robust(r,location,scale)
    cut=max(2.0_dp,abs(log(max(alpha,1.0e-12_dp))))
    out=min(location+cut*scale,max(location-cut*scale,r))
  end subroutine clean_boudt

  pure subroutine level_from_returns(r,start_value,levels)
    real(dp),intent(in)::r(:),start_value
    real(dp),intent(out)::levels(:)
    integer::i
    if(size(levels)==0) return
    levels(1)=start_value
    do i=1,min(size(r),size(levels)-1); levels(i+1)=levels(i)*(1.0_dp+r(i)); end do
  end subroutine level_from_returns

  pure subroutine convert_return_frequency(r,periods_per_output,out,nout,geometric)
    real(dp),intent(in)::r(:)
    integer,intent(in)::periods_per_output
    real(dp),intent(out)::out(:)
    integer,intent(out)::nout
    logical,intent(in),optional::geometric
    logical::g
    integer::i,a,b
    g=.true.; if(present(geometric)) g=geometric
    if(periods_per_output<=0) then; nout=0; return; end if
    nout=min(size(out),size(r)/periods_per_output)
    do i=1,nout
      a=(i-1)*periods_per_output+1; b=i*periods_per_output
      out(i)=cumulative_return(r(a:b),g)
    end do
  end subroutine convert_return_frequency

  subroutine portfolio_returns(r,weights,result,rebalance_every,initial_value,transaction_cost)
    real(dp),intent(in)::r(:,:),weights(:)
    type(portfolio_result),intent(out)::result
    integer,intent(in),optional::rebalance_every
    real(dp),intent(in),optional::initial_value,transaction_cost
    integer::n,p,t,k
    real(dp)::value,cost,pret,total,turn
    real(dp),allocatable::asset_values(:),target(:),bw(:),ew(:)
    n=size(r,1); p=size(r,2)
    allocate(result%returns(n),result%wealth(n),result%contributions(n,p))
    allocate(result%begin_weights(n,p),result%end_weights(n,p),result%turnover(n))
    allocate(asset_values(p),target(p),bw(p),ew(p))
    target=weights/max(sum(weights),tiny(1.0_dp)); value=1.0_dp
    if(present(initial_value)) value=initial_value
    cost=0.0_dp; if(present(transaction_cost)) cost=transaction_cost
    k=0; if(present(rebalance_every)) k=max(0,rebalance_every)
    asset_values=value*target
    do t=1,n
      total=sum(asset_values); if(total<=tiny(1.0_dp)) total=tiny(1.0_dp)
      bw=asset_values/total; result%begin_weights(t,:)=bw
      result%contributions(t,:)=bw*r(t,:); pret=sum(result%contributions(t,:))
      asset_values=asset_values*(1.0_dp+r(t,:)); value=sum(asset_values)
      result%returns(t)=pret; result%wealth(t)=value
      if(value>tiny(1.0_dp)) then; ew=asset_values/value; else; ew=target; end if
      result%end_weights(t,:)=ew; turn=0.0_dp
      if(k>0) then
        if(mod(t,k)==0 .and. t<n) then
          turn=0.5_dp*sum(abs(target-ew)); value=value*(1.0_dp-cost*turn)
          asset_values=value*target
        end if
      end if
      result%turnover(t)=turn
    end do
  end subroutine portfolio_returns
end module returns_mod
