! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of fExtremes.
! Copyright (C) 2026. This program is free software: you can redistribute it
! and/or modify it under the terms of the GNU General Public License as
! published by the Free Software Foundation, either version 2 of the License,
! or (at your option) any later version.
module fextremes_optimize
  use fextremes_kinds, only: dp
  implicit none
  private
  public :: objective_function, nelder_mead_bounded, numerical_hessian, invert_matrix

  abstract interface
    function objective_function(x) result(f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: f
    end function objective_function
  end interface
contains
  subroutine nelder_mead_bounded(fun, start, lower, upper, best, fbest, converged, iterations, evaluations, &
      max_iter, tolerance)
    procedure(objective_function) :: fun
    real(dp), intent(in) :: start(:), lower(:), upper(:)
    real(dp), intent(out) :: best(size(start)), fbest
    logical, intent(out) :: converged
    integer, intent(out) :: iterations, evaluations
    integer, intent(in), optional :: max_iter
    real(dp), intent(in), optional :: tolerance
    integer :: n, i, j, imax, isec, ilo, mit
    real(dp) :: tol, alpha, gamma, rho, sigma, fr, fe, fc, convspread
    real(dp), allocatable :: simplex(:,:), fv(:), centroid(:), xr(:), xe(:), xc(:)
    n = size(start); mit = 3000; tol = 1.0e-8_dp
    if (present(max_iter)) mit = max_iter
    if (present(tolerance)) tol = tolerance
    allocate(simplex(n,n+1),fv(n+1),centroid(n),xr(n),xe(n),xc(n))
    simplex(:,1) = min(max(start,lower),upper)
    do j=2,n+1
      simplex(:,j)=simplex(:,1)
      i=j-1
      simplex(i,j)=min(upper(i),max(lower(i),simplex(i,j)+0.05_dp*max(1.0_dp,abs(simplex(i,j)))))
      if (abs(simplex(i,j)-simplex(i,1)) <= epsilon(1.0_dp)) &
        simplex(i,j)=0.5_dp*(lower(i)+upper(i))
    end do
    do j=1,n+1; fv(j)=fun(simplex(:,j)); end do
    evaluations=n+1; converged=.false.; alpha=1.0_dp; gamma=2.0_dp; rho=0.5_dp; sigma=0.5_dp
    do iterations=1,mit
      ilo=minloc(fv,dim=1); imax=maxloc(fv,dim=1)
      isec=ilo
      do j=1,n+1
        if (j/=imax) then
          if (isec==imax .or. fv(j)>fv(isec)) isec=j
        end if
      end do
      convspread=maxval(abs(fv-fv(ilo)))/(1.0_dp+abs(fv(ilo)))
      if (convspread < tol .and. &
          maxval(abs(simplex-spread(simplex(:,ilo),2,n+1))) < &
          sqrt(tol)*(1.0_dp+maxval(abs(simplex(:,ilo))))) then
        converged=.true.; exit
      end if
      centroid=0.0_dp
      do j=1,n+1; if(j/=imax) centroid=centroid+simplex(:,j); end do
      centroid=centroid/real(n,dp)
      xr=min(max(centroid+alpha*(centroid-simplex(:,imax)),lower),upper); fr=fun(xr); evaluations=evaluations+1
      if (fr<fv(ilo)) then
        xe=min(max(centroid+gamma*(xr-centroid),lower),upper); fe=fun(xe); evaluations=evaluations+1
        if(fe<fr) then; simplex(:,imax)=xe; fv(imax)=fe; else; simplex(:,imax)=xr; fv(imax)=fr; end if
      else if (fr<fv(isec)) then
        simplex(:,imax)=xr; fv(imax)=fr
      else
        if(fr<fv(imax)) then; xc=min(max(centroid+rho*(xr-centroid),lower),upper)
        else; xc=min(max(centroid+rho*(simplex(:,imax)-centroid),lower),upper); end if
        fc=fun(xc); evaluations=evaluations+1
        if(fc<min(fr,fv(imax))) then
          simplex(:,imax)=xc; fv(imax)=fc
        else
          do j=1,n+1
            if(j/=ilo) then
              simplex(:,j)=min(max(simplex(:,ilo)+sigma*(simplex(:,j)-simplex(:,ilo)),lower),upper)
              fv(j)=fun(simplex(:,j)); evaluations=evaluations+1
            end if
          end do
        end if
      end if
    end do
    ilo=minloc(fv,dim=1); best=simplex(:,ilo); fbest=fv(ilo)
    if(iterations>mit) iterations=mit
  end subroutine nelder_mead_bounded

  subroutine numerical_hessian(fun, x, hess)
    procedure(objective_function) :: fun
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: hess(size(x),size(x))
    real(dp), allocatable :: xp(:), xm(:), xpp(:), xpm(:), xmp(:), xmm(:), step(:)
    real(dp) :: f0
    integer :: i,j,n
    n=size(x); allocate(xp(n),xm(n),xpp(n),xpm(n),xmp(n),xmm(n),step(n))
    step=epsilon(1.0_dp)**0.25_dp*max(1.0_dp,abs(x)); f0=fun(x)
    hess=0.0_dp
    do i=1,n
      xp=x; xm=x; xp(i)=xp(i)+step(i); xm(i)=xm(i)-step(i)
      hess(i,i)=(fun(xp)-2.0_dp*f0+fun(xm))/(step(i)**2)
      do j=i+1,n
        xpp=x; xpm=x; xmp=x; xmm=x
        xpp(i)=xpp(i)+step(i); xpp(j)=xpp(j)+step(j)
        xpm(i)=xpm(i)+step(i); xpm(j)=xpm(j)-step(j)
        xmp(i)=xmp(i)-step(i); xmp(j)=xmp(j)+step(j)
        xmm(i)=xmm(i)-step(i); xmm(j)=xmm(j)-step(j)
        hess(i,j)=(fun(xpp)-fun(xpm)-fun(xmp)+fun(xmm))/(4.0_dp*step(i)*step(j))
        hess(j,i)=hess(i,j)
      end do
    end do
  end subroutine numerical_hessian

  subroutine invert_matrix(a, ainv, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(size(a,1),size(a,2))
    logical, intent(out) :: ok
    real(dp), allocatable :: aug(:,:), rowtmp(:)
    real(dp) :: pivot, factor
    integer :: n,i,k,p
    n=size(a,1); ok=.false.; ainv=0.0_dp
    if(size(a,2)/=n) return
    allocate(aug(n,2*n),rowtmp(2*n)); aug(:,1:n)=a; aug(:,n+1:2*n)=0.0_dp
    do i=1,n; aug(i,n+i)=1.0_dp; end do
    do k=1,n
      p=k
      do i=k+1,n; if(abs(aug(i,k))>abs(aug(p,k))) p=i; end do
      if(abs(aug(p,k))<1.0e-12_dp) return
      if(p/=k) then; rowtmp=aug(k,:); aug(k,:)=aug(p,:); aug(p,:)=rowtmp; end if
      pivot=aug(k,k); aug(k,:)=aug(k,:)/pivot
      do i=1,n
        if(i/=k) then; factor=aug(i,k); aug(i,:)=aug(i,:)-factor*aug(k,:); end if
      end do
    end do
    ainv=aug(:,n+1:2*n); ok=.true.
  end subroutine invert_matrix
end module fextremes_optimize
