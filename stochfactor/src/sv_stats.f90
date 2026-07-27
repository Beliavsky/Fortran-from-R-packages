! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013-2026 the original stochvol and factorstochvol authors
! and the modern Fortran translation contributors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 2, or (at your option) any later version.
! This file is distributed without any warranty; see LICENSE for details.
module sv_stats
  use sv_kinds, only : dp, pi
  implicit none
  private
  public :: mean1, sd1, quantile_sorted, log_returns, normal_logpdf, normal_cdf
  public :: student_logpdf_std, safe_logsumexp, clamp
contains
  pure real(dp) function mean1(x) result(v)
    real(dp), intent(in) :: x(:); v=sum(x)/real(size(x),dp)
  end function mean1
  pure real(dp) function sd1(x) result(v)
    real(dp), intent(in) :: x(:); real(dp)::m
    m=mean1(x); v=sqrt(sum((x-m)**2)/real(max(1,size(x)-1),dp))
  end function sd1
  pure real(dp) function clamp(x,lo,hi) result(v)
    real(dp),intent(in)::x,lo,hi; v=min(hi,max(lo,x))
  end function clamp
  pure real(dp) function normal_logpdf(x,mu,sd) result(v)
    real(dp),intent(in)::x,mu,sd
    v=-0.5_dp*log(2.0_dp*pi)-log(sd)-0.5_dp*((x-mu)/sd)**2
  end function normal_logpdf
  pure real(dp) function normal_cdf(x) result(v)
    real(dp),intent(in)::x; v=0.5_dp*erfc(-x/sqrt(2.0_dp))
  end function normal_cdf
  pure real(dp) function student_logpdf_std(x,nu) result(v)
    real(dp),intent(in)::x,nu
    v=log_gamma(0.5_dp*(nu+1.0_dp))-log_gamma(0.5_dp*nu)-0.5_dp*log(pi*(nu-2.0_dp)) &
      -0.5_dp*(nu+1.0_dp)*log(1.0_dp+x*x/(nu-2.0_dp))
  end function student_logpdf_std
  pure real(dp) function safe_logsumexp(x) result(v)
    real(dp),intent(in)::x(:); real(dp)::m
    m=maxval(x); v=m+log(sum(exp(x-m)))
  end function safe_logsumexp
  subroutine log_returns(prices,returns,demean,standardize)
    real(dp),intent(in)::prices(:,:)
    real(dp),allocatable,intent(out)::returns(:,:)
    logical,intent(in),optional::demean,standardize
    logical::dm,st
    integer::i,j,n,m
    real(dp)::mu,s
    n=size(prices,1);m=size(prices,2);allocate(returns(n-1,m))
    do j=1,m; do i=1,n-1; returns(i,j)=log(prices(i+1,j)/prices(i,j)); end do; end do
    dm=.false.;st=.false.;if(present(demean))dm=demean;if(present(standardize))st=standardize
    do j=1,m
      if(dm) then; mu=mean1(returns(:,j)); returns(:,j)=returns(:,j)-mu; end if
      if(st) then; s=sd1(returns(:,j)); if(s>0.0_dp)returns(:,j)=returns(:,j)/s; end if
    end do
  end subroutine log_returns
  real(dp) function quantile_sorted(x,p) result(q)
    real(dp),intent(in)::x(:),p
    real(dp)::h,w
    integer::i,n
    n=size(x);h=1.0_dp+(n-1)*min(1.0_dp,max(0.0_dp,p));i=floor(h);w=h-i
    if(i>=n) then;q=x(n);else;q=(1.0_dp-w)*x(i)+w*x(i+1);end if
  end function quantile_sorted
end module sv_stats
