! SPDX-License-Identifier: GPL-2.0-or-later
module tmvtnorm_optimize
  use mvtnorm_kinds, only : dp
  implicit none
  private
  public :: nm_result_t, nelder_mead, scalar_objective

  type :: nm_result_t
    real(dp), allocatable :: par(:)
    real(dp) :: value = huge(1.0_dp)
    integer :: iterations = 0
    integer :: evaluations = 0
    integer :: convergence = 1
  end type nm_result_t

  abstract interface
    function scalar_objective(x) result(v)
      import :: dp
      real(dp), intent(in) :: x(:)
      real(dp) :: v
    end function scalar_objective
  end interface

contains

  function nelder_mead(fn, x0, maxit, reltol, step) result(out)
    procedure(scalar_objective) :: fn
    real(dp), intent(in) :: x0(:)
    integer, intent(in), optional :: maxit
    real(dp), intent(in), optional :: reltol, step
    type(nm_result_t) :: out
    integer :: n, miter, iter, j, ilo, ihi, inhi
    real(dp) :: tol, stp, fr, fe, fc, fbest, fspread, xspread
    real(dp), allocatable :: x(:,:), f(:), centroid(:), xr(:), xe(:), xc(:)
    real(dp), parameter :: alpha=1.0_dp, gamma=2.0_dp, rho=0.5_dp, shrink=0.5_dp

    n = size(x0)
    miter = 2000
    if (present(maxit)) miter = maxit
    tol = 1.0e-8_dp
    if (present(reltol)) tol = reltol
    stp = 0.05_dp
    if (present(step)) stp = step
    allocate(x(n,n+1), f(n+1), centroid(n), xr(n), xe(n), xc(n))
    x(:,1) = x0
    do j=1,n
      x(:,j+1)=x0
      x(j,j+1)=x(j,j+1)+stp*max(1.0_dp,abs(x0(j)))
    end do
    do j=1,n+1
      f(j)=safe_eval(fn,x(:,j))
    end do
    out%evaluations=n+1

    do iter=1,miter
      call extrema(f,ilo,ihi,inhi)
      fbest=f(ilo)
      fspread=maxval(abs(f-fbest))/max(1.0_dp,abs(fbest))
      xspread=0.0_dp
      do j=1,n+1
        xspread=max(xspread,maxval(abs(x(:,j)-x(:,ilo)))/max(1.0_dp,maxval(abs(x(:,ilo)))))
      end do
      if (fspread <= tol .and. xspread <= sqrt(tol)) then
        out%convergence=0
        exit
      end if

      centroid=0.0_dp
      do j=1,n+1
        if (j /= ihi) centroid=centroid+x(:,j)
      end do
      centroid=centroid/real(n,dp)

      xr=centroid+alpha*(centroid-x(:,ihi))
      fr=safe_eval(fn,xr)
      out%evaluations=out%evaluations+1
      if (fr < f(ilo)) then
        xe=centroid+gamma*(xr-centroid)
        fe=safe_eval(fn,xe)
        out%evaluations=out%evaluations+1
        if (fe < fr) then
          x(:,ihi)=xe
          f(ihi)=fe
        else
          x(:,ihi)=xr
          f(ihi)=fr
        end if
      else if (fr < f(inhi)) then
        x(:,ihi)=xr
        f(ihi)=fr
      else
        if (fr < f(ihi)) then
          xc=centroid+rho*(xr-centroid)
        else
          xc=centroid+rho*(x(:,ihi)-centroid)
        end if
        fc=safe_eval(fn,xc)
        out%evaluations=out%evaluations+1
        if (fc < min(fr,f(ihi))) then
          x(:,ihi)=xc
          f(ihi)=fc
        else
          do j=1,n+1
            if (j /= ilo) then
              x(:,j)=x(:,ilo)+shrink*(x(:,j)-x(:,ilo))
              f(j)=safe_eval(fn,x(:,j))
              out%evaluations=out%evaluations+1
            end if
          end do
        end if
      end if
    end do
    call extrema(f,ilo,ihi,inhi)
    out%par=x(:,ilo)
    out%value=f(ilo)
    out%iterations=min(iter,miter)
  end function nelder_mead

  real(dp) function safe_eval(fn,x) result(v)
    procedure(scalar_objective) :: fn
    real(dp), intent(in) :: x(:)
    v=fn(x)
    if (.not.(v < huge(1.0_dp))) v=huge(1.0_dp)/100.0_dp
  end function safe_eval

  subroutine extrema(f,ilo,ihi,inhi)
    real(dp), intent(in) :: f(:)
    integer, intent(out) :: ilo, ihi, inhi
    integer :: j
    ilo=minloc(f,dim=1)
    ihi=maxloc(f,dim=1)
    inhi=ilo
    do j=1,size(f)
      if (j==ihi) cycle
      if (inhi==ihi .or. f(j)>f(inhi)) inhi=j
    end do
  end subroutine extrema

end module tmvtnorm_optimize
