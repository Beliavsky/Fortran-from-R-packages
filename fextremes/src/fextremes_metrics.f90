! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of fExtremes.
! Copyright (C) 2026. This program is free software: you can redistribute it
! and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 2 of the License,
! or (at your option) any later version.
module fextremes_metrics
  use fextremes_kinds, only: dp
  use fextremes_stats, only: mean_value
  implicit none
  private
  public :: ema_filter, riskmetrics_volatility
contains
  subroutine ema_filter(x,lambda,ema,startup)
    real(dp),intent(in)::x(:),lambda
    real(dp),intent(out)::ema(size(x))
    integer,intent(in),optional::startup
    integer::i,nstart
    real(dp)::lam
    lam=lambda
    if(lam>=1.0_dp) lam=2.0_dp/(lam+1.0_dp)
    lam=max(0.0_dp,min(1.0_dp,lam))
    nstart=max(1,min(size(x),floor(2.0_dp/max(lam,epsilon(1.0_dp)))))
    if(present(startup)) then
      if(startup>0) nstart=max(1,min(size(x),startup))
    end if
    if(lam<=epsilon(1.0_dp)) then
      ema=mean_value(x)
      return
    end if
    ema(1)=mean_value(x(:nstart))
    do i=2,size(x); ema(i)=lam*x(i)+(1.0_dp-lam)*ema(i-1); end do
  end subroutine ema_filter

  subroutine riskmetrics_volatility(returns,decay,volatility,startup)
    real(dp),intent(in)::returns(:),decay
    real(dp),intent(out)::volatility(size(returns))
    integer,intent(in),optional::startup
    real(dp),allocatable::variance(:)
    allocate(variance(size(returns)))
    if(present(startup)) then
      call ema_filter(returns**2,1.0_dp-decay,variance,startup)
    else
      call ema_filter(returns**2,1.0_dp-decay,variance)
    end if
    volatility=sqrt(max(variance,0.0_dp))
  end subroutine riskmetrics_volatility
end module fextremes_metrics
