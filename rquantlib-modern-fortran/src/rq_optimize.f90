! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 OpenAI
! This file is part of a modern Fortran translation of RQuantLib.
! It is free software under GNU GPL version 2 or any later version.
module rq_optimize
  use rq_kinds, only: dp, huge_penalty
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none
  private
  public :: nelder_mead_bounded, optimize_result
  type :: optimize_result
    real(dp), allocatable :: x(:)
    real(dp) :: value = huge_penalty
    integer :: iterations = 0
    integer :: evaluations = 0
    integer :: status = 1
  end type optimize_result
  abstract interface
    function vector_objective(x) result(f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: f
    end function vector_objective
  end interface
contains
  subroutine nelder_mead_bounded(fun, start, lower, upper, result, tol, max_iter)
    procedure(vector_objective) :: fun
    real(dp), intent(in) :: start(:), lower(:), upper(:)
    type(optimize_result), intent(out) :: result
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: max_iter
    integer :: n, j, it, imax, imin, isec, nmax
    real(dp) :: eps, alpha, gamma, rho, sigma, fr, fe, fc, spread
    real(dp), allocatable :: p(:,:), f(:), centroid(:), xr(:), xe(:), xc(:)
    n=size(start); eps=1.0e-8_dp; nmax=2000
    if (present(tol)) eps=tol
    if (present(max_iter)) nmax=max_iter
    allocate(p(n,n+1),f(n+1),centroid(n),xr(n),xe(n),xc(n),result%x(n))
    p(:,1)=min(max(start,lower),upper)
    do j=1,n
      p(:,j+1)=p(:,1)
      p(j,j+1)=min(upper(j),max(lower(j),p(j,1)+0.05_dp*max(1.0_dp,abs(p(j,1)))))
      if (abs(p(j,j+1)-p(j,1)) <= epsilon(1.0_dp)) &
        p(j,j+1)=0.5_dp*(lower(j)+upper(j))
    end do
    do j=1,n+1
      f(j)=safe_eval(p(:,j)); result%evaluations=result%evaluations+1
    end do
    alpha=1.0_dp; gamma=2.0_dp; rho=0.5_dp; sigma=0.5_dp
    do it=1,nmax
      call extrema(f,imin,imax,isec)
      spread=maxval(abs(f-f(imin)))
      if (spread <= eps*(1.0_dp+abs(f(imin)))) exit
      centroid=0.0_dp
      do j=1,n+1
        if (j/=imax) centroid=centroid+p(:,j)
      end do
      centroid=centroid/real(n,dp)
      xr=min(max(centroid+alpha*(centroid-p(:,imax)),lower),upper)
      fr=safe_eval(xr); result%evaluations=result%evaluations+1
      if (fr < f(imin)) then
        xe=min(max(centroid+gamma*(xr-centroid),lower),upper)
        fe=safe_eval(xe); result%evaluations=result%evaluations+1
        if (fe<fr) then; p(:,imax)=xe; f(imax)=fe; else; p(:,imax)=xr; f(imax)=fr; end if
      else if (fr < f(isec)) then
        p(:,imax)=xr; f(imax)=fr
      else
        if (fr < f(imax)) then
          xc=min(max(centroid+rho*(xr-centroid),lower),upper)
        else
          xc=min(max(centroid-rho*(centroid-p(:,imax)),lower),upper)
        end if
        fc=safe_eval(xc); result%evaluations=result%evaluations+1
        if (fc < min(fr,f(imax))) then
          p(:,imax)=xc; f(imax)=fc
        else
          do j=1,n+1
            if (j==imin) cycle
            p(:,j)=min(max(p(:,imin)+sigma*(p(:,j)-p(:,imin)),lower),upper)
            f(j)=safe_eval(p(:,j)); result%evaluations=result%evaluations+1
          end do
        end if
      end if
    end do
    call extrema(f,imin,imax,isec)
    result%x=p(:,imin); result%value=f(imin); result%iterations=it
    if (it<=nmax) result%status=0
  contains
    function safe_eval(x) result(v)
      real(dp), intent(in) :: x(:)
      real(dp) :: v
      v=fun(x)
      if (ieee_is_nan(v) .or. abs(v)>=huge_penalty) v=huge_penalty
    end function safe_eval
    subroutine extrema(v,ilo,ihi,isecond)
      real(dp), intent(in) :: v(:)
      integer, intent(out) :: ilo,ihi,isecond
      integer :: k
      ilo=1; ihi=1
      do k=2,size(v)
        if(v(k)<v(ilo)) ilo=k
        if(v(k)>v(ihi)) ihi=k
      end do
      isecond=ilo
      do k=1,size(v)
        if(k==ihi) cycle
        if(isecond==ihi .or. v(k)>v(isecond)) isecond=k
      end do
    end subroutine extrema
  end subroutine nelder_mead_bounded
end module rq_optimize
