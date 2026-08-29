! SPDX-License-Identifier: GPL-2.0-or-later
module nlme_optimize
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use nlme_kinds, only : dp
  use nlme_status, only : NLME_SUCCESS, NLME_INVALID_ARGUMENT, NLME_MAX_ITER, NLME_NONFINITE
  implicit none
  private
  public :: objective_function, nelder_mead, finite_difference_gradient, finite_difference_hessian

  abstract interface
    function objective_function(x, context) result(value)
      import dp
      real(dp), intent(in) :: x(:)
      class(*), intent(inout) :: context
      real(dp) :: value
    end function objective_function
  end interface
contains
  subroutine nelder_mead(fun, context, x, value, status, iterations, max_iter, tolerance, initial_step, verbose)
    procedure(objective_function) :: fun
    class(*), intent(inout) :: context
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out) :: value
    integer, intent(out) :: status, iterations
    integer, intent(in), optional :: max_iter
    real(dp), intent(in), optional :: tolerance, initial_step
    logical, intent(in), optional :: verbose
    integer :: n, miter, i, j, ilo, ihi, inhi
    real(dp) :: tol, step, fr, fe, fc, spread_f, spread_x
    real(dp), allocatable :: simplex(:,:), f(:), centroid(:), xr(:), xe(:), xc(:), temp(:)
    logical :: trace
    real(dp), parameter :: alpha=1.0_dp, gamma=2.0_dp, rho=0.5_dp, sigma=0.5_dp

    n=size(x)
    miter=1000
    if (present(max_iter)) miter=max_iter
    tol=1.0e-8_dp
    if (present(tolerance)) tol=max(tolerance,epsilon(1.0_dp))
    step=0.1_dp
    if (present(initial_step)) step=max(initial_step,sqrt(epsilon(1.0_dp)))
    trace=.false.
    if (present(verbose)) trace=verbose
    if (n==0 .or. miter<1) then
      value=huge(1.0_dp)
      status=NLME_INVALID_ARGUMENT
      iterations=0
      return
    end if
    allocate(simplex(n,n+1),f(n+1),centroid(n),xr(n),xe(n),xc(n),temp(n))
    simplex(:,1)=x
    do j=2,n+1
      simplex(:,j)=x
      i=j-1
      simplex(i,j)=x(i)+step*max(1.0_dp,abs(x(i)))
    end do
    do j=1,n+1
      f(j)=fun(simplex(:,j),context)
      if (.not.ieee_is_finite(f(j))) f(j)=huge(1.0_dp)/100.0_dp
    end do
    status=NLME_MAX_ITER
    do iterations=1,miter
      ilo=minloc(f,dim=1)
      ihi=maxloc(f,dim=1)
      inhi=ilo
      do j=1,n+1
        if (j==ihi) cycle
        if (inhi==ihi .or. f(j)>f(inhi)) inhi=j
      end do
      spread_f=maxval(abs(f-f(ilo)))/max(1.0_dp,abs(f(ilo)))
      spread_x=0.0_dp
      do j=1,n+1
        spread_x=max(spread_x,maxval(abs(simplex(:,j)-simplex(:,ilo))))
      end do
      spread_x=spread_x/max(1.0_dp,maxval(abs(simplex(:,ilo))))
      if (spread_f<=tol .and. spread_x<=sqrt(tol)) then
        status=NLME_SUCCESS
        exit
      end if
      centroid=(sum(simplex,dim=2)-simplex(:,ihi))/real(n,dp)
      xr=centroid+alpha*(centroid-simplex(:,ihi))
      fr=fun(xr,context)
      if (.not.ieee_is_finite(fr)) fr=huge(1.0_dp)/100.0_dp
      if (fr<f(ilo)) then
        xe=centroid+gamma*(xr-centroid)
        fe=fun(xe,context)
        if (.not.ieee_is_finite(fe)) fe=huge(1.0_dp)/100.0_dp
        if (fe<fr) then
          simplex(:,ihi)=xe
          f(ihi)=fe
        else
          simplex(:,ihi)=xr
          f(ihi)=fr
        end if
      else if (fr<f(inhi)) then
        simplex(:,ihi)=xr
        f(ihi)=fr
      else
        if (fr<f(ihi)) then
          xc=centroid+rho*(xr-centroid)
        else
          xc=centroid+rho*(simplex(:,ihi)-centroid)
        end if
        fc=fun(xc,context)
        if (.not.ieee_is_finite(fc)) fc=huge(1.0_dp)/100.0_dp
        if (fc<min(fr,f(ihi))) then
          simplex(:,ihi)=xc
          f(ihi)=fc
        else
          temp=simplex(:,ilo)
          do j=1,n+1
            if (j==ilo) cycle
            simplex(:,j)=temp+sigma*(simplex(:,j)-temp)
            f(j)=fun(simplex(:,j),context)
            if (.not.ieee_is_finite(f(j))) f(j)=huge(1.0_dp)/100.0_dp
          end do
        end if
      end if
      if (trace .and. mod(iterations,50)==0) then
        write(*,'(a,i0,a,es14.6)') 'nelder_mead iteration ',iterations,': ',minval(f)
      end if
    end do
    ilo=minloc(f,dim=1)
    x=simplex(:,ilo)
    value=f(ilo)
    if (.not.ieee_is_finite(value)) status=NLME_NONFINITE
  end subroutine nelder_mead

  subroutine finite_difference_gradient(fun,context,x,gradient,status,step)
    procedure(objective_function) :: fun
    class(*), intent(inout) :: context
    real(dp), intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: gradient(:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: step
    real(dp), allocatable :: xp(:),xm(:)
    real(dp) :: h,fp,fm
    integer :: i
    h=epsilon(1.0_dp)**(1.0_dp/3.0_dp)
    if (present(step)) h=step
    allocate(gradient(size(x)),xp(size(x)),xm(size(x)))
    do i=1,size(x)
      xp=x
      xm=x
      xp(i)=xp(i)+h*max(1.0_dp,abs(x(i)))
      xm(i)=xm(i)-h*max(1.0_dp,abs(x(i)))
      fp=fun(xp,context)
      fm=fun(xm,context)
      if (.not.ieee_is_finite(fp) .or. .not.ieee_is_finite(fm)) then
        gradient=0.0_dp
        status=NLME_NONFINITE
        return
      end if
      gradient(i)=(fp-fm)/(xp(i)-xm(i))
    end do
    status=NLME_SUCCESS
  end subroutine finite_difference_gradient

  subroutine finite_difference_hessian(fun,context,x,hessian,status,step)
    procedure(objective_function) :: fun
    class(*), intent(inout) :: context
    real(dp), intent(in) :: x(:)
    real(dp), allocatable, intent(out) :: hessian(:,:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: step
    real(dp), allocatable :: xpp(:),xpm(:),xmp(:),xmm(:)
    real(dp) :: h,hi,hj,fpp,fpm,fmp,fmm,f0,fp,fm
    integer :: i,j,n
    h=epsilon(1.0_dp)**0.25_dp
    if (present(step)) h=step
    n=size(x)
    allocate(hessian(n,n),xpp(n),xpm(n),xmp(n),xmm(n))
    hessian=0.0_dp
    f0=fun(x,context)
    if (.not.ieee_is_finite(f0)) then
    status=NLME_NONFINITE
    return
    end if
    do i=1,n
      hi=h*max(1.0_dp,abs(x(i)))
      xpp=x
      xmm=x
      xpp(i)=x(i)+hi
      xmm(i)=x(i)-hi
      fp=fun(xpp,context)
      fm=fun(xmm,context)
      if (.not.ieee_is_finite(fp) .or. .not.ieee_is_finite(fm)) then
      status=NLME_NONFINITE
      return
      end if
      hessian(i,i)=(fp-2.0_dp*f0+fm)/(hi*hi)
      do j=1,i-1
        hj=h*max(1.0_dp,abs(x(j)))
        xpp=x
        xpm=x
        xmp=x
        xmm=x
        xpp(i)=x(i)+hi
        xpp(j)=x(j)+hj
        xpm(i)=x(i)+hi
        xpm(j)=x(j)-hj
        xmp(i)=x(i)-hi
        xmp(j)=x(j)+hj
        xmm(i)=x(i)-hi
        xmm(j)=x(j)-hj
        fpp=fun(xpp,context)
        fpm=fun(xpm,context)
        fmp=fun(xmp,context)
        fmm=fun(xmm,context)
        if (any(.not.ieee_is_finite([fpp,fpm,fmp,fmm]))) then
        status=NLME_NONFINITE
        return
        end if
        hessian(i,j)=(fpp-fpm-fmp+fmm)/(4.0_dp*hi*hj)
        hessian(j,i)=hessian(i,j)
      end do
    end do
    status=NLME_SUCCESS
  end subroutine finite_difference_hessian
end module nlme_optimize
