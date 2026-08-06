! SPDX-License-Identifier: GPL-2.0-only
module tsd_math
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use ghyp_kinds, only : dp, pi, sqrt_two_pi
  use ghyp_special, only : normal_pdf, normal_cdf, normal_quantile, regularized_beta, student_cdf
  implicit none
  private

  public :: normal_pdf, normal_cdf, normal_quantile, student_cdf
  public :: student_pdf, student_quantile, regularized_gamma_p, gamma_quantile
  public :: signum, heaviside, safe_log, clamp_probability
  public :: mean_value, sample_variance, sample_sd, quantile_type7, sort_real
  public :: invert_matrix, solve_linear, finite_vector, log_sum_exp_vector
  public :: linear_interpolate, trapz, lower_bound_index

contains

  pure real(dp) function signum(x) result(value)
    real(dp), intent(in) :: x
    if (x > 0.0_dp) then
      value = 1.0_dp
    else if (x < 0.0_dp) then
      value = -1.0_dp
    else
      value = 0.0_dp
    end if
  end function signum

  pure real(dp) function heaviside(x, at_zero) result(value)
    real(dp), intent(in) :: x
    real(dp), intent(in), optional :: at_zero
    real(dp) :: a
    a = 0.5_dp
    if (present(at_zero)) a = at_zero
    if (x > 0.0_dp) then
      value = 1.0_dp
    else if (x < 0.0_dp) then
      value = 0.0_dp
    else
      value = a
    end if
  end function heaviside

  pure real(dp) function safe_log(x) result(value)
    real(dp), intent(in) :: x
    value = log(max(x, tiny(1.0_dp)))
  end function safe_log

  pure real(dp) function clamp_probability(p) result(value)
    real(dp), intent(in) :: p
    value = min(1.0_dp, max(0.0_dp, p))
  end function clamp_probability

  pure real(dp) function student_pdf(x, nu) result(value)
    real(dp), intent(in) :: x, nu
    if (nu <= 0.0_dp) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      value = exp(log_gamma(0.5_dp*(nu+1.0_dp)) - log_gamma(0.5_dp*nu) - &
                  0.5_dp*log(pi*nu) - 0.5_dp*(nu+1.0_dp)*log(1.0_dp+x*x/nu))
    end if
  end function student_pdf

  real(dp) function student_quantile(p, nu) result(value)
    real(dp), intent(in) :: p, nu
    real(dp) :: lo, hi, mid, fmid
    integer :: iter
    if (p <= 0.0_dp) then
      value = -huge(1.0_dp)
      return
    else if (p >= 1.0_dp) then
      value = huge(1.0_dp)
      return
    else if (nu <= 0.0_dp) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    lo = -1.0_dp
    hi = 1.0_dp
    do while (student_cdf(lo,nu) > p)
      lo = 2.0_dp*lo
      if (lo < -1.0e8_dp) exit
    end do
    do while (student_cdf(hi,nu) < p)
      hi = 2.0_dp*hi
      if (hi > 1.0e8_dp) exit
    end do
    do iter = 1, 160
      mid = 0.5_dp*(lo+hi)
      fmid = student_cdf(mid,nu)
      if (fmid < p) then
        lo = mid
      else
        hi = mid
      end if
      if (abs(hi-lo) <= 2.0e-12_dp*max(1.0_dp,abs(mid))) exit
    end do
    value = 0.5_dp*(lo+hi)
  end function student_quantile

  pure real(dp) function regularized_gamma_p(a, x) result(value)
    real(dp), intent(in) :: a, x
    real(dp), parameter :: eps = 3.0e-14_dp, fpmin = 1.0e-300_dp
    real(dp) :: ap, del, sumv, b, c, d, h, an
    integer :: n
    if (a <= 0.0_dp .or. x < 0.0_dp) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    else if (abs(x) <= tiny(1.0_dp)) then
      value = 0.0_dp
      return
    end if
    if (x < a + 1.0_dp) then
      ap = a
      sumv = 1.0_dp/a
      del = sumv
      do n = 1, 400
        ap = ap + 1.0_dp
        del = del*x/ap
        sumv = sumv + del
        if (abs(del) <= abs(sumv)*eps) exit
      end do
      value = sumv*exp(-x+a*log(x)-log_gamma(a))
    else
      b = x + 1.0_dp - a
      c = 1.0_dp/fpmin
      d = 1.0_dp/max(abs(b),fpmin)
      if (b < 0.0_dp) d = -d
      h = d
      do n = 1, 400
        an = -real(n,dp)*(real(n,dp)-a)
        b = b + 2.0_dp
        d = an*d+b
        if (abs(d) < fpmin) d = fpmin
        c = b+an/c
        if (abs(c) < fpmin) c = fpmin
        d = 1.0_dp/d
        del = d*c
        h = h*del
        if (abs(del-1.0_dp) <= eps) exit
      end do
      value = 1.0_dp-exp(-x+a*log(x)-log_gamma(a))*h
    end if
    value = clamp_probability(value)
  end function regularized_gamma_p

  real(dp) function gamma_quantile(p, shape, scale) result(value)
    real(dp), intent(in) :: p, shape
    real(dp), intent(in), optional :: scale
    real(dp) :: sc, lo, hi, mid, pmid
    integer :: iter
    sc = 1.0_dp
    if (present(scale)) sc = scale
    if (p <= 0.0_dp) then
      value = 0.0_dp
      return
    else if (p >= 1.0_dp) then
      value = huge(1.0_dp)
      return
    else if (shape <= 0.0_dp .or. sc <= 0.0_dp) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    lo = 0.0_dp
    hi = max(shape + 8.0_dp*sqrt(shape) + 8.0_dp, 1.0_dp)
    do while (regularized_gamma_p(shape,hi) < p)
      hi = 2.0_dp*hi
      if (hi > 1.0e12_dp) exit
    end do
    do iter = 1, 180
      mid = 0.5_dp*(lo+hi)
      pmid = regularized_gamma_p(shape,mid)
      if (pmid < p) then
        lo = mid
      else
        hi = mid
      end if
      if (abs(hi-lo) <= 2.0e-12_dp*max(1.0_dp,abs(mid))) exit
    end do
    value = sc*0.5_dp*(lo+hi)
  end function gamma_quantile

  pure real(dp) function mean_value(x) result(value)
    real(dp), intent(in) :: x(:)
    if (size(x) == 0) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      value = sum(x)/real(size(x),dp)
    end if
  end function mean_value

  pure real(dp) function sample_variance(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    if (size(x) < 2) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      m = mean_value(x)
      value = sum((x-m)**2)/real(size(x)-1,dp)
    end if
  end function sample_variance

  pure real(dp) function sample_sd(x) result(value)
    real(dp), intent(in) :: x(:)
    value = sqrt(max(sample_variance(x),0.0_dp))
  end function sample_sd

  subroutine sort_real(x)
    real(dp), intent(inout) :: x(:)
    real(dp) :: key
    integer :: i, j
    do i = 2, size(x)
      key = x(i)
      j = i-1
      do while (j >= 1)
        if (x(j) <= key) exit
        x(j+1) = x(j)
        j = j-1
      end do
      x(j+1) = key
    end do
  end subroutine sort_real

  real(dp) function quantile_type7(x, p) result(value)
    real(dp), intent(in) :: x(:), p
    real(dp), allocatable :: y(:)
    real(dp) :: h, frac
    integer :: j, n
    n = size(x)
    if (n == 0 .or. p < 0.0_dp .or. p > 1.0_dp) then
      value = ieee_value(0.0_dp,ieee_quiet_nan)
      return
    end if
    allocate(y(n)); y = x; call sort_real(y)
    if (p <= 0.0_dp) then
      value = y(1)
    else if (p >= 1.0_dp) then
      value = y(n)
    else
      h = 1.0_dp + real(n-1,dp)*p
      j = int(floor(h))
      frac = h-real(j,dp)
      if (j >= n) then
        value = y(n)
      else
        value = (1.0_dp-frac)*y(j)+frac*y(j+1)
      end if
    end if
  end function quantile_type7

  subroutine solve_linear(a, b, x, ok)
    real(dp), intent(in) :: a(:, :), b(:)
    real(dp), allocatable, intent(out) :: x(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: aug(:, :), row(:)
    real(dp) :: pivot, factor
    integer :: n, i, k, p
    n = size(b)
    allocate(x(n))
    if (size(a,1) /= n .or. size(a,2) /= n) then
      x = 0.0_dp; ok = .false.; return
    end if
    allocate(aug(n,n+1),row(n+1)); aug(:,1:n)=a; aug(:,n+1)=b
    ok = .true.
    do k=1,n
      p=k
      do i=k+1,n
        if(abs(aug(i,k))>abs(aug(p,k))) p=i
      end do
      if(abs(aug(p,k)) <= 1.0e-14_dp*max(1.0_dp,maxval(abs(a)))) then
        ok=.false.; x=0.0_dp; return
      end if
      if(p/=k) then
        row=aug(k,:); aug(k,:)=aug(p,:); aug(p,:)=row
      end if
      pivot=aug(k,k); aug(k,:)=aug(k,:)/pivot
      do i=1,n
        if(i==k) cycle
        factor=aug(i,k)
        aug(i,:)=aug(i,:)-factor*aug(k,:)
      end do
    end do
    x=aug(:,n+1)
  end subroutine solve_linear

  subroutine invert_matrix(a, inverse, ok)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: inverse(:, :)
    logical, intent(out) :: ok
    real(dp), allocatable :: aug(:, :), row(:)
    real(dp) :: pivot, factor
    integer :: n, i, k, p
    n=size(a,1); allocate(inverse(n,n))
    if(size(a,2)/=n) then
      inverse=0.0_dp; ok=.false.; return
    end if
    allocate(aug(n,2*n),row(2*n)); aug=0.0_dp; aug(:,1:n)=a
    do i=1,n; aug(i,n+i)=1.0_dp; end do
    ok=.true.
    do k=1,n
      p=k
      do i=k+1,n
        if(abs(aug(i,k))>abs(aug(p,k))) p=i
      end do
      if(abs(aug(p,k))<=1.0e-13_dp*max(1.0_dp,maxval(abs(a)))) then
        inverse=0.0_dp;ok=.false.;return
      end if
      if(p/=k) then; row=aug(k,:);aug(k,:)=aug(p,:);aug(p,:)=row;end if
      pivot=aug(k,k);aug(k,:)=aug(k,:)/pivot
      do i=1,n
        if(i==k) cycle
        factor=aug(i,k);aug(i,:)=aug(i,:)-factor*aug(k,:)
      end do
    end do
    inverse=aug(:,n+1:2*n)
  end subroutine invert_matrix

  pure logical function finite_vector(x) result(ok)
    real(dp), intent(in) :: x(:)
    ok = all(ieee_is_finite(x))
  end function finite_vector

  pure real(dp) function log_sum_exp_vector(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    if (size(x)==0) then
      value=-huge(1.0_dp)
    else
      m=maxval(x)
      value=m+log(sum(exp(x-m)))
    end if
  end function log_sum_exp_vector

  pure real(dp) function linear_interpolate(x, grid_x, grid_y) result(value)
    real(dp), intent(in) :: x, grid_x(:), grid_y(:)
    integer :: j
    real(dp) :: t
    if(size(grid_x)/=size(grid_y) .or. size(grid_x)==0) then
      value=ieee_value(0.0_dp,ieee_quiet_nan);return
    end if
    if(x<=grid_x(1)) then;value=grid_y(1);return;end if
    if(x>=grid_x(size(grid_x))) then;value=grid_y(size(grid_y));return;end if
    j=lower_bound_index(grid_x,x)
    t=(x-grid_x(j))/(grid_x(j+1)-grid_x(j))
    value=(1.0_dp-t)*grid_y(j)+t*grid_y(j+1)
  end function linear_interpolate

  pure integer function lower_bound_index(grid, x) result(index)
    real(dp), intent(in) :: grid(:), x
    integer :: lo, hi, mid
    lo=1;hi=size(grid)-1
    do while(lo<hi)
      mid=(lo+hi+1)/2
      if(grid(mid)<=x) then;lo=mid;else;hi=mid-1;end if
    end do
    index=max(1,min(size(grid)-1,lo))
  end function lower_bound_index

  pure real(dp) function trapz(x, y) result(value)
    real(dp), intent(in) :: x(:), y(:)
    integer :: i
    value=0.0_dp
    if(size(x)/=size(y) .or. size(x)<2) return
    do i=1,size(x)-1
      value=value+0.5_dp*(x(i+1)-x(i))*(y(i+1)+y(i))
    end do
  end function trapz

end module tsd_math
