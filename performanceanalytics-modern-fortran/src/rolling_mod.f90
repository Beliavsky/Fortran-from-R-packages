! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of PerformanceAnalytics.
! Copyright (C) original PerformanceAnalytics authors and translation contributors.
! Distributed under GNU GPL version 2 or, at your option, any later version.
module rolling_mod
  use kinds_mod, only: dp
  use statistics_mod, only: mean_value, sd_value, quantile_type7, correlation_value
  use returns_mod, only: cumulative_return, annualized_return
  use risk_mod, only: value_at_risk, expected_shortfall
  implicit none
  private
  public :: rolling_statistic, expanding_statistic, rolling_correlation
contains
  subroutine rolling_statistic(r,width,statistic,out,nout,scale,p,method)
    real(dp),intent(in)::r(:)
    integer,intent(in)::width
    character(len=*),intent(in)::statistic
    real(dp),intent(out)::out(:)
    integer,intent(out)::nout
    real(dp),intent(in),optional::scale,p
    character(len=*),intent(in),optional::method
    real(dp)::sc,prob
    character(len=24)::meth
    integer::i
    sc=252.0_dp;if(present(scale))sc=scale
    prob=0.95_dp;if(present(p))prob=p
    meth='modified';if(present(method))meth=method
    if(width<=0 .or. size(r)<width)then;nout=0;return;end if
    nout=min(size(out),size(r)-width+1)
    do i=1,nout
      select case(trim(adjustl(statistic)))
      case('mean');out(i)=mean_value(r(i:i+width-1))
      case('sd','stddev');out(i)=sd_value(r(i:i+width-1))
      case('cumulative');out(i)=cumulative_return(r(i:i+width-1))
      case('annualized');out(i)=annualized_return(r(i:i+width-1),sc)
      case('var');out(i)=value_at_risk(r(i:i+width-1),prob,meth)
      case('es');out(i)=expected_shortfall(r(i:i+width-1),prob,meth)
      case('median');out(i)=quantile_type7(r(i:i+width-1),0.5_dp)
      case default;out(i)=mean_value(r(i:i+width-1))
      end select
    end do
  end subroutine rolling_statistic

  subroutine expanding_statistic(r,min_obs,statistic,out,nout,scale,p,method)
    real(dp),intent(in)::r(:)
    integer,intent(in)::min_obs
    character(len=*),intent(in)::statistic
    real(dp),intent(out)::out(:)
    integer,intent(out)::nout
    real(dp),intent(in),optional::scale,p
    character(len=*),intent(in),optional::method
    real(dp)::sc,prob
    character(len=24)::meth
    integer::i
    sc=252.0_dp;if(present(scale))sc=scale
    prob=0.95_dp;if(present(p))prob=p
    meth='modified';if(present(method))meth=method
    if(min_obs<=0 .or. size(r)<min_obs)then;nout=0;return;end if
    nout=min(size(out),size(r)-min_obs+1)
    do i=1,nout
      select case(trim(adjustl(statistic)))
      case('mean');out(i)=mean_value(r(:min_obs+i-1))
      case('sd','stddev');out(i)=sd_value(r(:min_obs+i-1))
      case('cumulative');out(i)=cumulative_return(r(:min_obs+i-1))
      case('annualized');out(i)=annualized_return(r(:min_obs+i-1),sc)
      case('var');out(i)=value_at_risk(r(:min_obs+i-1),prob,meth)
      case('es');out(i)=expected_shortfall(r(:min_obs+i-1),prob,meth)
      case default;out(i)=mean_value(r(:min_obs+i-1))
      end select
    end do
  end subroutine expanding_statistic

  subroutine rolling_correlation(x,y,width,out,nout)
    real(dp),intent(in)::x(:),y(:)
    integer,intent(in)::width
    real(dp),intent(out)::out(:)
    integer,intent(out)::nout
    integer::i,n
    n=min(size(x),size(y))
    if(width<=0 .or. n<width)then;nout=0;return;end if
    nout=min(size(out),n-width+1)
    do i=1,nout;out(i)=correlation_value(x(i:i+width-1),y(i:i+width-1));end do
  end subroutine rolling_correlation
end module rolling_mod
