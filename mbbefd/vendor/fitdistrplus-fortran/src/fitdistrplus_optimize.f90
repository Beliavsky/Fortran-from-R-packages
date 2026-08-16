! SPDX-License-Identifier: GPL-2.0-or-later
module fitdistrplus_optimize
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fitdistrplus_kinds, only : dp
  use fitdistrplus_types, only : fit_success, fit_invalid_argument, fit_no_convergence
  implicit none
  private

  public :: nelder_mead, parameters_to_unconstrained, unconstrained_to_parameters
  public :: transformation_jacobian

  abstract interface
    function objective_function(x, context) result(value)
      import dp
      real(dp), intent(in) :: x(:)
      class(*), intent(inout) :: context
      real(dp) :: value
    end function objective_function
  end interface

contains

  subroutine nelder_mead(fun, context, x, fval, status, iterations, max_iterations, tolerance, scale)
    procedure(objective_function) :: fun
    class(*), intent(inout) :: context
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out) :: fval
    integer, intent(out) :: status, iterations
    integer, intent(in), optional :: max_iterations
    real(dp), intent(in), optional :: tolerance, scale
    real(dp), allocatable :: simplex(:, :), values(:), centroid(:), xr(:), xe(:), xc(:)
    real(dp) :: alpha, gamma, rho, sigma, tol, step, fr, fe, fc
    real(dp) :: spread_x, spread_f
    real(dp), allocatable :: tmpx(:)
    integer :: n, j, worst, second_worst, best, maxit

    n = size(x)
    if (n == 0) then
      fval = fun(x,context); status = fit_invalid_argument; iterations = 0; return
    end if
    maxit = 2000
    if (present(max_iterations)) maxit = max_iterations
    tol = 1.0e-8_dp
    if (present(tolerance)) tol = tolerance
    step = 0.10_dp
    if (present(scale)) step = scale
    alpha = 1.0_dp; gamma = 2.0_dp; rho = 0.5_dp; sigma = 0.5_dp

    allocate(simplex(n,n+1),values(n+1),centroid(n),xr(n),xe(n),xc(n),tmpx(n))
    simplex(:,1) = x
    do j = 1, n
      simplex(:,j+1) = x
      simplex(j,j+1) = x(j) + step*max(abs(x(j)),1.0_dp)
    end do
    do j = 1, n+1
      values(j) = safe_objective(fun,context,simplex(:,j))
    end do

    status = fit_no_convergence
    do iterations = 1, maxit
      call sort_simplex(simplex,values)
      best = 1
      second_worst = n
      worst = n+1
      spread_f = maxval(abs(values-values(best)))
      spread_x = 0.0_dp
      do j = 2, n+1
        spread_x = max(spread_x,maxval(abs(simplex(:,j)-simplex(:,best))))
      end do
      if (spread_f <= tol*(1.0_dp+abs(values(best))) .and. &
          spread_x <= sqrt(tol)*(1.0_dp+maxval(abs(simplex(:,best))))) then
        status = fit_success
        exit
      end if

      centroid = sum(simplex(:,1:n),dim=2)/real(n,dp)
      xr = centroid + alpha*(centroid-simplex(:,worst))
      fr = safe_objective(fun,context,xr)
      if (fr < values(best)) then
        xe = centroid + gamma*(xr-centroid)
        fe = safe_objective(fun,context,xe)
        if (fe < fr) then
          simplex(:,worst)=xe; values(worst)=fe
        else
          simplex(:,worst)=xr; values(worst)=fr
        end if
      else if (fr < values(second_worst)) then
        simplex(:,worst)=xr; values(worst)=fr
      else
        if (fr < values(worst)) then
          xc = centroid + rho*(xr-centroid)
        else
          xc = centroid - rho*(centroid-simplex(:,worst))
        end if
        fc = safe_objective(fun,context,xc)
        if (fc < min(fr,values(worst))) then
          simplex(:,worst)=xc; values(worst)=fc
        else
          do j = 2, n+1
            simplex(:,j)=simplex(:,best)+sigma*(simplex(:,j)-simplex(:,best))
            values(j)=safe_objective(fun,context,simplex(:,j))
          end do
        end if
      end if
    end do
    call sort_simplex(simplex,values)
    x = simplex(:,1)
    fval = values(1)
    if (iterations > maxit) iterations = maxit
  end subroutine nelder_mead

  function safe_objective(fun,context,x) result(value)
    procedure(objective_function) :: fun
    class(*), intent(inout) :: context
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    value=fun(x,context)
    if (.not.ieee_is_finite(value)) value=huge(1.0_dp)/100.0_dp
  end function safe_objective

  subroutine sort_simplex(simplex,values)
    real(dp), intent(inout) :: simplex(:, :), values(:)
    real(dp), allocatable :: column(:)
    real(dp) :: key
    integer :: i, j
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

  subroutine parameters_to_unconstrained(par,lower,upper,z,status)
    real(dp), intent(in) :: par(:),lower(:),upper(:)
    real(dp), allocatable, intent(out) :: z(:)
    integer, intent(out) :: status
    integer :: i,n
    real(dp) :: p
    n=size(par)
    allocate(z(n))
    if(size(lower)/=n .or. size(upper)/=n .or. any(lower>=upper)) then
      z=0.0_dp; status=fit_invalid_argument; return
    end if
    do i=1,n
      p=par(i)
      if(p<=lower(i) .or. p>=upper(i)) then
        if(is_finite_bound(lower(i)) .and. is_finite_bound(upper(i))) then
          p=min(max(p,lower(i)+1.0e-8_dp*(upper(i)-lower(i))), &
            upper(i)-1.0e-8_dp*(upper(i)-lower(i)))
        else if(is_finite_bound(lower(i))) then
          p=max(p,lower(i)+1.0e-8_dp*max(1.0_dp,abs(lower(i))))
        else if(is_finite_bound(upper(i))) then
          p=min(p,upper(i)-1.0e-8_dp*max(1.0_dp,abs(upper(i))))
        end if
      end if
      if(is_finite_bound(lower(i)) .and. is_finite_bound(upper(i))) then
        z(i)=log((p-lower(i))/(upper(i)-p))
      else if(is_finite_bound(lower(i))) then
        z(i)=log(p-lower(i))
      else if(is_finite_bound(upper(i))) then
        z(i)=log(upper(i)-p)
      else
        z(i)=p
      end if
    end do
    status=fit_success
  end subroutine parameters_to_unconstrained

  subroutine unconstrained_to_parameters(z,lower,upper,par)
    real(dp), intent(in) :: z(:),lower(:),upper(:)
    real(dp), intent(out) :: par(:)
    integer :: i
    real(dp) :: s
    do i=1,size(z)
      if(is_finite_bound(lower(i)) .and. is_finite_bound(upper(i))) then
        if(z(i)>=0.0_dp) then
          s=1.0_dp/(1.0_dp+exp(-min(z(i),700.0_dp)))
        else
          s=exp(max(z(i),-700.0_dp))/(1.0_dp+exp(max(z(i),-700.0_dp)))
        end if
        par(i)=lower(i)+(upper(i)-lower(i))*s
      else if(is_finite_bound(lower(i))) then
        par(i)=lower(i)+exp(min(z(i),700.0_dp))
      else if(is_finite_bound(upper(i))) then
        par(i)=upper(i)-exp(min(z(i),700.0_dp))
      else
        par(i)=z(i)
      end if
    end do
  end subroutine unconstrained_to_parameters

  subroutine transformation_jacobian(z,lower,upper,jacobian)
    real(dp), intent(in) :: z(:),lower(:),upper(:)
    real(dp), allocatable, intent(out) :: jacobian(:, :)
    real(dp) :: s
    integer :: i,n
    n=size(z); allocate(jacobian(n,n)); jacobian=0.0_dp
    do i=1,n
      if(is_finite_bound(lower(i)) .and. is_finite_bound(upper(i))) then
        if(z(i)>=0.0_dp) then
          s=1.0_dp/(1.0_dp+exp(-min(z(i),700.0_dp)))
        else
          s=exp(max(z(i),-700.0_dp))/(1.0_dp+exp(max(z(i),-700.0_dp)))
        end if
        jacobian(i,i)=(upper(i)-lower(i))*s*(1.0_dp-s)
      else if(is_finite_bound(lower(i))) then
        jacobian(i,i)=exp(min(z(i),700.0_dp))
      else if(is_finite_bound(upper(i))) then
        jacobian(i,i)=-exp(min(z(i),700.0_dp))
      else
        jacobian(i,i)=1.0_dp
      end if
    end do
  end subroutine transformation_jacobian

  pure logical function is_finite_bound(value)
    real(dp), intent(in) :: value
    is_finite_bound = abs(value) < 0.25_dp*huge(1.0_dp)
  end function is_finite_bound

end module fitdistrplus_optimize
