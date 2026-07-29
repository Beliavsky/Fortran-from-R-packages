! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2016 Marius Hofert, Kurt Hornik and Alexander J. McNeil
module qrmtools_optimization
  use qrmtools_kinds, only : dp
  implicit none
  private
  public :: nelder_mead
  abstract interface
    function objective_function(x) result(value)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: value
    end function objective_function
  end interface
contains
  subroutine nelder_mead(f,start,optimum,best_value,iterations,evaluations,converged,&
      max_iterations,tolerance)
    procedure(objective_function) :: f
    real(dp), intent(in) :: start(:)
    real(dp), allocatable, intent(out) :: optimum(:)
    real(dp), intent(out) :: best_value
    integer, intent(out) :: iterations,evaluations
    logical, intent(out) :: converged
    integer, intent(in), optional :: max_iterations
    real(dp), intent(in), optional :: tolerance
    real(dp), allocatable :: simplex(:,:),values(:),centroid(:),xr(:),xe(:),xc(:)
    real(dp) :: step,fr,fe,fc,tol
    integer :: n,i,j,maxit
    n=size(start); maxit=2000; if(present(max_iterations))maxit=max_iterations
    tol=1.0e-9_dp; if(present(tolerance))tol=tolerance
    allocate(simplex(n,n+1),values(n+1),centroid(n),xr(n),xe(n),xc(n))
    simplex(:,1)=start
    do i=1,n
      simplex(:,i+1)=start
      step=0.05_dp*max(abs(start(i)),1.0_dp)
      simplex(i,i+1)=simplex(i,i+1)+step
    end do
    do i=1,n+1; values(i)=f(simplex(:,i)); end do
    evaluations=n+1; converged=.false.
    do iterations=1,maxit
      call sort_simplex(simplex,values)
      if(maxval(abs(simplex(:,2:n+1)-spread(simplex(:,1),2,n)))/&
          max(1.0_dp,maxval(abs(simplex(:,1))))<tol .and. &
          maxval(abs(values(2:n+1)-values(1)))<tol) then
        converged=.true.; exit
      end if
      centroid=sum(simplex(:,1:n),dim=2)/real(n,dp)
      xr=centroid+(centroid-simplex(:,n+1)); fr=f(xr); evaluations=evaluations+1
      if(fr<values(1)) then
        xe=centroid+2.0_dp*(xr-centroid); fe=f(xe); evaluations=evaluations+1
        if(fe<fr) then; simplex(:,n+1)=xe; values(n+1)=fe
        else; simplex(:,n+1)=xr; values(n+1)=fr; end if
      else if(fr<values(n)) then
        simplex(:,n+1)=xr; values(n+1)=fr
      else
        if(fr<values(n+1)) then; xc=centroid+0.5_dp*(xr-centroid)
        else; xc=centroid+0.5_dp*(simplex(:,n+1)-centroid); end if
        fc=f(xc); evaluations=evaluations+1
        if(fc<min(fr,values(n+1))) then
          simplex(:,n+1)=xc; values(n+1)=fc
        else
          do j=2,n+1
            simplex(:,j)=simplex(:,1)+0.5_dp*(simplex(:,j)-simplex(:,1))
            values(j)=f(simplex(:,j))
          end do
          evaluations=evaluations+n
        end if
      end if
    end do
    call sort_simplex(simplex,values); optimum=simplex(:,1); best_value=values(1)
  end subroutine nelder_mead

  subroutine sort_simplex(simplex,values)
    real(dp), intent(inout) :: simplex(:,:),values(:)
    integer :: i,j
    real(dp) :: key
    real(dp), allocatable :: column(:)
    allocate(column(size(simplex,1)))
    do i=2,size(values)
      key=values(i); column=simplex(:,i); j=i-1
      do while(j>=1)
        if(values(j)<=key) exit
        values(j+1)=values(j); simplex(:,j+1)=simplex(:,j); j=j-1
      end do
      values(j+1)=key; simplex(:,j+1)=column
    end do
  end subroutine sort_simplex
end module qrmtools_optimization
