! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from highfrequency 1.0.2 by Kris Boudt, Jonathan Cornelissen,
! Scott Payseur, Onno Kleen, Emil Sjoerup, and contributors.
module highfrequency_optimize
  use highfrequency_kinds, only: dp
  implicit none
  private
  public :: nelder_mead_bounded

  abstract interface
    function objective_function(x) result(value)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: value
    end function objective_function
  end interface

contains

  subroutine nelder_mead_bounded(objective, start, lower, upper, solution, value, converged, max_iter, tolerance)
    procedure(objective_function) :: objective
    real(dp), intent(in) :: start(:), lower(:), upper(:)
    real(dp), intent(out) :: solution(size(start)), value
    logical, intent(out) :: converged
    integer, intent(in), optional :: max_iter
    real(dp), intent(in), optional :: tolerance
    real(dp), allocatable :: simplex(:,:), f(:), centroid(:), xr(:), xe(:), xc(:)
    real(dp) :: alpha, gamma, rho, sigma, tol, spread, step, fr, fe, fc
    integer :: n, iter, maxit, i, j, best, worst, second_worst, tmp
    integer, allocatable :: order(:)

    n=size(start)
    allocate(simplex(n,n+1),f(n+1),centroid(n),xr(n),xe(n),xc(n),order(n+1))
    alpha=1.0_dp
    gamma=2.0_dp
    rho=0.5_dp
    sigma=0.5_dp
    maxit=2000
    if(present(max_iter)) maxit=max_iter
    tol=1.0e-9_dp
    if(present(tolerance)) tol=tolerance
    simplex(:,1)=min(upper,max(lower,start))
    do j=1,n
      simplex(:,j+1)=simplex(:,1)
      step=0.05_dp*max(1.0_dp,abs(start(j)))
      if(upper(j)>lower(j)) step=max(step,0.05_dp*(upper(j)-lower(j)))
      simplex(j,j+1)=min(upper(j),max(lower(j),simplex(j,j+1)+step))
      if(simplex(j,j+1)==simplex(j,1)) simplex(j,j+1)=min(upper(j),max(lower(j),simplex(j,j+1)-2.0_dp*step))
    end do
    do j=1,n+1
      f(j)=objective(simplex(:,j))
    end do
    converged=.false.
    do iter=1,maxit
      order=[(j,j=1,n+1)]
      do i=2,n+1
        tmp=order(i)
        j=i-1
        do while(j>=1)
          if(f(order(j))<=f(tmp)) exit
          order(j+1)=order(j)
          j=j-1
        end do
        order(j+1)=tmp
      end do
      best=order(1)
      worst=order(n+1)
      second_worst=order(n)
      spread=maxval(abs(simplex-spread_ref(simplex(:,best),n+1)))
      if(spread<=tol .and. maxval(abs(f-f(best)))<=sqrt(tol))then
        converged=.true.
        exit
      end if
      centroid=0.0_dp
      do i=1,n
        centroid=centroid+simplex(:,order(i))
      end do
      centroid=centroid/real(n,dp)
      xr=min(upper,max(lower,centroid+alpha*(centroid-simplex(:,worst))))
      fr=objective(xr)
      if(f(best)<=fr .and. fr<f(second_worst))then
        simplex(:,worst)=xr
        f(worst)=fr
      else if(fr<f(best))then
        xe=min(upper,max(lower,centroid+gamma*(xr-centroid)))
        fe=objective(xe)
        if(fe<fr)then
          simplex(:,worst)=xe
          f(worst)=fe
        else
          simplex(:,worst)=xr
          f(worst)=fr
        end if
      else
        xc=min(upper,max(lower,centroid+rho*(simplex(:,worst)-centroid)))
        fc=objective(xc)
        if(fc<f(worst))then
          simplex(:,worst)=xc
          f(worst)=fc
        else
          do i=2,n+1
            simplex(:,order(i))=min(upper,max(lower,simplex(:,best)+sigma*(simplex(:,order(i))-simplex(:,best))))
            f(order(i))=objective(simplex(:,order(i)))
          end do
        end if
      end if
    end do
    best=minloc(f,dim=1)
    solution=simplex(:,best)
    value=f(best)
  contains
    pure function spread_ref(x,m) result(a)
      real(dp),intent(in)::x(:)
      integer,intent(in)::m
      real(dp)::a(size(x),m)
      integer::k
      do k=1,m
        a(:,k)=x
      end do
    end function spread_ref
  end subroutine nelder_mead_bounded

end module highfrequency_optimize
