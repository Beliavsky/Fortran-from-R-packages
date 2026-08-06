! SPDX-License-Identifier: GPL-3.0-only
module tscopula_math
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use tscopula_kinds, only : dp, pi, log_two_pi
  use tscopula_status, only : tsc_error, tsc_success, tsc_invalid_input, &
    tsc_numerical_failure, tsc_not_converged, clear_error, set_error
  implicit none
  private

  abstract interface
    function scalar_function(x, context) result(value)
      import dp
      real(dp), intent(in) :: x
      class(*), intent(inout) :: context
      real(dp) :: value
    end function scalar_function
    function vector_objective(x, context) result(value)
      import dp
      real(dp), intent(in) :: x(:)
      class(*), intent(inout) :: context
      real(dp) :: value
    end function vector_objective
  end interface

  type, public :: optimizer_result
    real(dp), allocatable :: par(:)
    real(dp) :: value = huge(1.0_dp)
    integer :: iterations = 0
    integer :: evaluations = 0
    integer :: convergence = 1
  end type optimizer_result

  public :: normal_pdf, normal_cdf, normal_quantile
  public :: student_pdf, student_cdf, student_quantile
  public :: beta_pdf, regularized_beta, beta_quantile
  public :: uniform_random, normal_random, student_random, set_seed
  public :: integrate_simpson, bisection_root, minimize_nelder_mead
  public :: solve_linear, invert_matrix, determinant_logabs, symmetrize
  public :: finite_hessian, safe_standard_errors, empirical_quantile_type7
  public :: log_sum_exp, clamp_probability, finite_all

contains

  elemental real(dp) function clamp_probability(p) result(q)
    real(dp), intent(in) :: p
    real(dp), parameter :: epsp = 2.0_dp*epsilon(1.0_dp)
    q = min(max(p, epsp), 1.0_dp-epsp)
  end function clamp_probability

  pure logical function finite_all(x) result(ok)
    real(dp), intent(in) :: x(:)
    integer :: i
    ok = .true.
    do i = 1, size(x)
      if (.not. ieee_is_finite(x(i))) then
        ok = .false.
        return
      end if
    end do
  end function finite_all

  elemental real(dp) function normal_pdf(x) result(value)
    real(dp), intent(in) :: x
    value = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
  end function normal_pdf

  elemental real(dp) function normal_cdf(x) result(value)
    real(dp), intent(in) :: x
    value = 0.5_dp*erfc(-x/sqrt(2.0_dp))
  end function normal_cdf

  elemental real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp), parameter :: a1 = -3.969683028665376e1_dp
    real(dp), parameter :: a2 =  2.209460984245205e2_dp
    real(dp), parameter :: a3 = -2.759285104469687e2_dp
    real(dp), parameter :: a4 =  1.383577518672690e2_dp
    real(dp), parameter :: a5 = -3.066479806614716e1_dp
    real(dp), parameter :: a6 =  2.506628277459239_dp
    real(dp), parameter :: b1 = -5.447609879822406e1_dp
    real(dp), parameter :: b2 =  1.615858368580409e2_dp
    real(dp), parameter :: b3 = -1.556989798598866e2_dp
    real(dp), parameter :: b4 =  6.680131188771972e1_dp
    real(dp), parameter :: b5 = -1.328068155288572e1_dp
    real(dp), parameter :: c1 = -7.784894002430293e-3_dp
    real(dp), parameter :: c2 = -3.223964580411365e-1_dp
    real(dp), parameter :: c3 = -2.400758277161838_dp
    real(dp), parameter :: c4 = -2.549732539343734_dp
    real(dp), parameter :: c5 =  4.374664141464968_dp
    real(dp), parameter :: c6 =  2.938163982698783_dp
    real(dp), parameter :: d1 =  7.784695709041462e-3_dp
    real(dp), parameter :: d2 =  3.224671290700398e-1_dp
    real(dp), parameter :: d3 =  2.445134137142996_dp
    real(dp), parameter :: d4 =  3.754408661907416_dp
    real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp-plow
    real(dp) :: q, r
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (p < plow) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
        ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    else if (p <= phigh) then
      q = p-0.5_dp
      r = q*q
      x = (((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6)*q / &
        (((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
    else
      q = sqrt(-2.0_dp*log(1.0_dp-p))
      x = -(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6) / &
        ((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
    end if
    if (p > 0.0_dp .and. p < 1.0_dp) then
      x = x - (normal_cdf(x)-p)/max(normal_pdf(x), tiny(1.0_dp))
    end if
  end function normal_quantile

  elemental real(dp) function beta_pdf(x, a, b) result(value)
    real(dp), intent(in) :: x, a, b
    if (a <= 0.0_dp .or. b <= 0.0_dp .or. x < 0.0_dp .or. x > 1.0_dp) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
    else if (abs(x) <= tiny(1.0_dp)) then
      if (a < 1.0_dp) then
        value = huge(1.0_dp)
      else if (abs(a-1.0_dp) <= epsilon(1.0_dp)) then
        value = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b))
      else
        value = 0.0_dp
      end if
    else if (abs(x-1.0_dp) <= tiny(1.0_dp)) then
      if (b < 1.0_dp) then
        value = huge(1.0_dp)
      else if (abs(b-1.0_dp) <= epsilon(1.0_dp)) then
        value = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b))
      else
        value = 0.0_dp
      end if
    else
      value = exp((a-1.0_dp)*log(x)+(b-1.0_dp)*log(1.0_dp-x) + &
        log_gamma(a+b)-log_gamma(a)-log_gamma(b))
    end if
  end function beta_pdf

  pure real(dp) function beta_continued_fraction(a, b, x) result(cf)
    real(dp), intent(in) :: a, b, x
    integer, parameter :: max_iter = 300
    real(dp), parameter :: fpmin = tiny(1.0_dp)/epsilon(1.0_dp)
    real(dp) :: aa, c, d, del, h, qab, qam, qap
    integer :: m, m2
    qab = a+b
    qap = a+1.0_dp
    qam = a-1.0_dp
    c = 1.0_dp
    d = 1.0_dp-qab*x/qap
    if (abs(d) < fpmin) d = fpmin
    d = 1.0_dp/d
    h = d
    do m = 1, max_iter
      m2 = 2*m
      aa = real(m,dp)*(b-real(m,dp))*x / &
        ((qam+real(m2,dp))*(a+real(m2,dp)))
      d = 1.0_dp+aa*d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp+aa/c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp/d
      h = h*d*c
      aa = -(a+real(m,dp))*(qab+real(m,dp))*x / &
        ((a+real(m2,dp))*(qap+real(m2,dp)))
      d = 1.0_dp+aa*d
      if (abs(d) < fpmin) d = fpmin
      c = 1.0_dp+aa/c
      if (abs(c) < fpmin) c = fpmin
      d = 1.0_dp/d
      del = d*c
      h = h*del
      if (abs(del-1.0_dp) <= 8.0_dp*epsilon(1.0_dp)) exit
    end do
    cf = h
  end function beta_continued_fraction

  elemental real(dp) function regularized_beta(x, a, b) result(value)
    real(dp), intent(in) :: x, a, b
    real(dp) :: bt
    if (a <= 0.0_dp .or. b <= 0.0_dp) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    if (x <= 0.0_dp) then
      value = 0.0_dp
      return
    else if (x >= 1.0_dp) then
      value = 1.0_dp
      return
    end if
    bt = exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b) + &
      a*log(x)+b*log(1.0_dp-x))
    if (x < (a+1.0_dp)/(a+b+2.0_dp)) then
      value = bt*beta_continued_fraction(a,b,x)/a
    else
      value = 1.0_dp-bt*beta_continued_fraction(b,a,1.0_dp-x)/b
    end if
    value = min(max(value,0.0_dp),1.0_dp)
  end function regularized_beta

  elemental real(dp) function beta_quantile(p, a, b) result(x)
    real(dp), intent(in) :: p, a, b
    real(dp) :: lo, hi, mid
    integer :: iter
    if (p <= 0.0_dp) then
      x = 0.0_dp
      return
    else if (p >= 1.0_dp) then
      x = 1.0_dp
      return
    else if (a <= 0.0_dp .or. b <= 0.0_dp) then
      x = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    lo = 0.0_dp
    hi = 1.0_dp
    do iter = 1, 100
      mid = 0.5_dp*(lo+hi)
      if (regularized_beta(mid,a,b) < p) then
        lo = mid
      else
        hi = mid
      end if
    end do
    x = 0.5_dp*(lo+hi)
  end function beta_quantile

  elemental real(dp) function student_pdf(x, nu) result(value)
    real(dp), intent(in) :: x, nu
    if (nu <= 0.0_dp) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
    else
      value = exp(log_gamma(0.5_dp*(nu+1.0_dp))-log_gamma(0.5_dp*nu) - &
        0.5_dp*log(nu*pi)-0.5_dp*(nu+1.0_dp)*log(1.0_dp+x*x/nu))
    end if
  end function student_pdf

  elemental real(dp) function student_cdf(x, nu) result(value)
    real(dp), intent(in) :: x, nu
    real(dp) :: ib
    if (nu <= 0.0_dp) then
      value = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    if (abs(x) <= tiny(1.0_dp)) then
      value = 0.5_dp
      return
    end if
    ib = regularized_beta(nu/(nu+x*x), 0.5_dp*nu, 0.5_dp)
    if (x > 0.0_dp) then
      value = 1.0_dp-0.5_dp*ib
    else
      value = 0.5_dp*ib
    end if
  end function student_cdf

  elemental real(dp) function student_quantile(p, nu) result(x)
    real(dp), intent(in) :: p, nu
    real(dp) :: lo, hi, mid
    integer :: iter
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
      return
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
      return
    else if (nu <= 0.0_dp) then
      x = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    lo = -1.0_dp
    hi = 1.0_dp
    do while (student_cdf(lo,nu) > p)
      lo = 2.0_dp*lo
    end do
    do while (student_cdf(hi,nu) < p)
      hi = 2.0_dp*hi
    end do
    do iter = 1, 100
      mid = 0.5_dp*(lo+hi)
      if (student_cdf(mid,nu) < p) then
        lo = mid
      else
        hi = mid
      end if
    end do
    x = 0.5_dp*(lo+hi)
  end function student_quantile

  subroutine set_seed(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: values(:)
    call random_seed(size=n)
    allocate(values(n))
    do i = 1, n
      values(i) = modulo(seed+104729*i, huge(1)-1)
      if (values(i) == 0) values(i) = i
    end do
    call random_seed(put=values)
  end subroutine set_seed

  real(dp) function uniform_random() result(value)
    call random_number(value)
    value = clamp_probability(value)
  end function uniform_random

  real(dp) function normal_random() result(value)
    real(dp) :: u1, u2
    u1 = uniform_random()
    u2 = uniform_random()
    value = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
  end function normal_random

  recursive real(dp) function gamma_random(shape) result(value)
    real(dp), intent(in) :: shape
    real(dp) :: d, c, x, v, u
    if (shape < 1.0_dp) then
      value = gamma_random(shape+1.0_dp)*uniform_random()**(1.0_dp/shape)
      return
    end if
    d = shape-1.0_dp/3.0_dp
    c = 1.0_dp/sqrt(9.0_dp*d)
    do
      do
        x = normal_random()
        v = 1.0_dp+c*x
        if (v > 0.0_dp) exit
      end do
      v = v*v*v
      u = uniform_random()
      if (u < 1.0_dp-0.0331_dp*x**4) exit
      if (log(u) < 0.5_dp*x*x+d*(1.0_dp-v+log(v))) exit
    end do
    value = d*v
  end function gamma_random

  real(dp) function student_random(nu) result(value)
    real(dp), intent(in) :: nu
    value = normal_random()/sqrt(2.0_dp*gamma_random(0.5_dp*nu)/nu)
  end function student_random

  recursive function integrate_simpson(func, a, b, context, tol, max_depth) result(value)
    procedure(scalar_function) :: func
    real(dp), intent(in) :: a, b
    class(*), intent(inout) :: context
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: max_depth
    real(dp) :: value, tolerance, fa, fb, fc, whole
    integer :: depth
    tolerance = 1.0e-9_dp
    if (present(tol)) tolerance = tol
    depth = 18
    if (present(max_depth)) depth = max_depth
    fa = func(a,context)
    fb = func(b,context)
    fc = func(0.5_dp*(a+b),context)
    whole = (b-a)*(fa+4.0_dp*fc+fb)/6.0_dp
    value = adaptive_simpson(func,a,b,fa,fb,fc,whole,tolerance,depth,context)
  end function integrate_simpson

  recursive function adaptive_simpson(func,a,b,fa,fb,fc,whole,tol,depth,context) result(value)
    procedure(scalar_function) :: func
    real(dp), intent(in) :: a,b,fa,fb,fc,whole,tol
    integer, intent(in) :: depth
    class(*), intent(inout) :: context
    real(dp) :: value, c, d, e, fd, fe, left, right
    c = 0.5_dp*(a+b)
    d = 0.5_dp*(a+c)
    e = 0.5_dp*(c+b)
    fd = func(d,context)
    fe = func(e,context)
    left = (c-a)*(fa+4.0_dp*fd+fc)/6.0_dp
    right = (b-c)*(fc+4.0_dp*fe+fb)/6.0_dp
    if (depth <= 0 .or. abs(left+right-whole) <= 15.0_dp*tol) then
      value = left+right+(left+right-whole)/15.0_dp
    else
      value = adaptive_simpson(func,a,c,fa,fc,fd,left,0.5_dp*tol,depth-1,context) + &
        adaptive_simpson(func,c,b,fc,fb,fe,right,0.5_dp*tol,depth-1,context)
    end if
  end function adaptive_simpson

  subroutine bisection_root(func, lower, upper, context, root, error, tol)
    procedure(scalar_function) :: func
    real(dp), intent(in) :: lower, upper
    class(*), intent(inout) :: context
    real(dp), intent(out) :: root
    type(tsc_error), intent(out) :: error
    real(dp), intent(in), optional :: tol
    real(dp) :: a,b,c,fa,fb,fc,tolerance
    integer :: iter
    call clear_error(error)
    tolerance = 1.0e-10_dp
    if (present(tol)) tolerance = tol
    a = lower
    b = upper
    fa = func(a,context)
    fb = func(b,context)
    if (.not. ieee_is_finite(fa) .or. .not. ieee_is_finite(fb) .or. fa*fb > 0.0_dp) then
      root = ieee_value(0.0_dp,ieee_quiet_nan)
      call set_error(error,tsc_invalid_input,'root is not bracketed')
      return
    end if
    do iter = 1, 200
      c = 0.5_dp*(a+b)
      fc = func(c,context)
      if (abs(fc) <= tolerance .or. abs(b-a) <= tolerance*(1.0_dp+abs(c))) exit
      if (fa*fc <= 0.0_dp) then
        b = c
        fb = fc
      else
        a = c
        fa = fc
      end if
    end do
    root = 0.5_dp*(a+b)
  end subroutine bisection_root

  subroutine solve_linear(a, b, x, error)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), allocatable, intent(out) :: x(:)
    type(tsc_error), intent(out) :: error
    real(dp), allocatable :: m(:,:), rhs(:), row(:)
    real(dp) :: pivot, factor, tmp
    integer :: i,k,n,p
    call clear_error(error)
    n = size(b)
    if (size(a,1) /= n .or. size(a,2) /= n) then
      allocate(x(0))
      call set_error(error,tsc_invalid_input,'linear system has incompatible dimensions')
      return
    end if
    allocate(m(n,n),rhs(n),row(n),x(n))
    m = a
    rhs = b
    do k = 1,n
      p = k-1+maxloc(abs(m(k:n,k)),dim=1)
      pivot = m(p,k)
      if (abs(pivot) <= 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(m)))) then
        x = ieee_value(0.0_dp,ieee_quiet_nan)
        call set_error(error,tsc_numerical_failure,'singular linear system')
        return
      end if
      if (p /= k) then
        row = m(k,:); m(k,:) = m(p,:); m(p,:) = row
        tmp = rhs(k); rhs(k) = rhs(p); rhs(p) = tmp
      end if
      do i = k+1,n
        factor = m(i,k)/m(k,k)
        m(i,k:n) = m(i,k:n)-factor*m(k,k:n)
        rhs(i) = rhs(i)-factor*rhs(k)
      end do
    end do
    do i = n,1,-1
      x(i) = (rhs(i)-dot_product(m(i,i+1:n),x(i+1:n)))/m(i,i)
    end do
  end subroutine solve_linear

  subroutine invert_matrix(a, ainv, error)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    type(tsc_error), intent(out) :: error
    real(dp), allocatable :: rhs(:), col(:)
    type(tsc_error) :: local_error
    integer :: i,n
    call clear_error(error)
    n = size(a,1)
    if (size(a,2) /= n) then
      allocate(ainv(0,0))
      call set_error(error,tsc_invalid_input,'matrix must be square')
      return
    end if
    allocate(ainv(n,n),rhs(n))
    do i = 1,n
      rhs = 0.0_dp
      rhs(i) = 1.0_dp
      call solve_linear(a,rhs,col,local_error)
      if (.not. local_error%ok()) then
        deallocate(ainv)
        allocate(ainv(0,0))
        error = local_error
        return
      end if
      ainv(:,i) = col
    end do
  end subroutine invert_matrix

  subroutine determinant_logabs(a, logabs, sign_det, error)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: logabs, sign_det
    type(tsc_error), intent(out) :: error
    real(dp), allocatable :: m(:,:), row(:)
    real(dp) :: pivot, factor
    integer :: i,k,n,p
    call clear_error(error)
    n = size(a,1)
    if (size(a,2) /= n) then
      call set_error(error,tsc_invalid_input,'matrix must be square')
      logabs = -huge(1.0_dp); sign_det = 0.0_dp
      return
    end if
    allocate(m(n,n),row(n)); m=a
    logabs=0.0_dp; sign_det=1.0_dp
    do k=1,n
      p=k-1+maxloc(abs(m(k:n,k)),dim=1)
      pivot=m(p,k)
      if (abs(pivot) <= tiny(1.0_dp)) then
        logabs=-huge(1.0_dp); sign_det=0.0_dp
        call set_error(error,tsc_numerical_failure,'singular matrix')
        return
      end if
      if (p/=k) then
        row=m(k,:);m(k,:)=m(p,:);m(p,:)=row;sign_det=-sign_det
      end if
      if (m(k,k)<0.0_dp) sign_det=-sign_det
      logabs=logabs+log(abs(m(k,k)))
      do i=k+1,n
        factor=m(i,k)/m(k,k)
        m(i,k+1:n)=m(i,k+1:n)-factor*m(k,k+1:n)
      end do
    end do
  end subroutine determinant_logabs

  subroutine symmetrize(a)
    real(dp), intent(inout) :: a(:,:)
    a = 0.5_dp*(a+transpose(a))
  end subroutine symmetrize

  subroutine minimize_nelder_mead(objective, start, lower, upper, context, result, &
      max_iter, tolerance)
    procedure(vector_objective) :: objective
    real(dp), intent(in) :: start(:), lower(:), upper(:)
    class(*), intent(inout) :: context
    type(optimizer_result), intent(out) :: result
    integer, intent(in), optional :: max_iter
    real(dp), intent(in), optional :: tolerance
    real(dp), allocatable :: simplex(:,:), f(:), centroid(:), xr(:), xe(:), xc(:)
    real(dp), allocatable :: temp_col(:)
    real(dp) :: alpha,gamma,rho,sigma,tol,simplex_spread,step
    integer :: j,n,it,maxit,best,worst,second_worst
    n=size(start)
    if (size(lower)/=n .or. size(upper)/=n) then
      allocate(result%par(0));result%value=huge(1.0_dp);result%convergence=2
      return
    end if
    maxit=2000;if(present(max_iter))maxit=max_iter
    tol=1.0e-8_dp;if(present(tolerance))tol=tolerance
    alpha=1.0_dp;gamma=2.0_dp;rho=0.5_dp;sigma=0.5_dp
    allocate(simplex(n,n+1),f(n+1),centroid(n),xr(n),xe(n),xc(n),temp_col(n))
    simplex(:,1)=min(max(start,lower),upper)
    do j=2,n+1
      simplex(:,j)=simplex(:,1)
      step=0.05_dp*max(1.0_dp,abs(simplex(j-1,1)))
      if (upper(j-1)>lower(j-1)) step=min(step,0.25_dp*(upper(j-1)-lower(j-1)))
      simplex(j-1,j)=min(max(simplex(j-1,j)+step,lower(j-1)),upper(j-1))
      if (abs(simplex(j-1,j)-simplex(j-1,1)) <= tiny(1.0_dp)) &
        simplex(j-1,j)=min(max(simplex(j-1,j)-step,lower(j-1)),upper(j-1))
    end do
    do j=1,n+1
      f(j)=objective(simplex(:,j),context)
      if (.not.ieee_is_finite(f(j))) f(j)=huge(1.0_dp)/100.0_dp
    end do
    result%evaluations=n+1
    do it=1,maxit
      call sort_simplex(simplex,f)
      best=1;worst=n+1;second_worst=n
      simplex_spread=maxval(abs(simplex-spread(simplex(:,best),2,n+1)))
      if (simplex_spread <= tol*(1.0_dp+maxval(abs(simplex(:,best)))) .and. &
          abs(f(worst)-f(best)) <= sqrt(tol)*(1.0_dp+abs(f(best)))) exit
      centroid=sum(simplex(:,1:n),dim=2)/real(n,dp)
      xr=min(max(centroid+alpha*(centroid-simplex(:,worst)),lower),upper)
      result%evaluations=result%evaluations+1
      if (objective(xr,context)<f(best)) then
        f(worst)=objective(xr,context)
        result%evaluations=result%evaluations+1
        xe=min(max(centroid+gamma*(xr-centroid),lower),upper)
        if (objective(xe,context)<f(worst)) then
          f(worst)=objective(xe,context);simplex(:,worst)=xe
        else
          simplex(:,worst)=xr
        end if
        result%evaluations=result%evaluations+2
      else
        f(worst)=objective(xr,context)
        result%evaluations=result%evaluations+1
        if (f(worst)<f(second_worst)) then
          simplex(:,worst)=xr
        else
          if (f(worst)<objective(simplex(:,worst),context)) then
            xc=min(max(centroid+rho*(xr-centroid),lower),upper)
          else
            xc=min(max(centroid-rho*(centroid-simplex(:,worst)),lower),upper)
          end if
          result%evaluations=result%evaluations+1
          if (objective(xc,context)<min(f(worst),objective(simplex(:,worst),context))) then
            simplex(:,worst)=xc
            f(worst)=objective(xc,context)
            result%evaluations=result%evaluations+3
          else
            do j=2,n+1
              simplex(:,j)=min(max(simplex(:,best)+sigma*(simplex(:,j)-simplex(:,best)),lower),upper)
              f(j)=objective(simplex(:,j),context)
              result%evaluations=result%evaluations+1
            end do
          end if
        end if
      end if
    end do
    call sort_simplex(simplex,f)
    allocate(result%par(n));result%par=simplex(:,1);result%value=f(1)
    result%iterations=min(it,maxit)
    result%convergence=merge(0,1,it<=maxit)
  contains
    subroutine sort_simplex(x,fx)
      real(dp),intent(inout)::x(:,:),fx(:)
      real(dp)::tf
      integer::ii,jj,kk
      do ii=1,size(fx)-1
        kk=ii
        do jj=ii+1,size(fx)
          if(fx(jj)<fx(kk))kk=jj
        end do
        if(kk/=ii)then
          tf=fx(ii);fx(ii)=fx(kk);fx(kk)=tf
          temp_col=x(:,ii);x(:,ii)=x(:,kk);x(:,kk)=temp_col
        end if
      end do
    end subroutine sort_simplex
  end subroutine minimize_nelder_mead

  subroutine finite_hessian(objective, x, context, hessian)
    procedure(vector_objective) :: objective
    real(dp), intent(in) :: x(:)
    class(*), intent(inout) :: context
    real(dp), allocatable, intent(out) :: hessian(:,:)
    real(dp), allocatable :: xp(:),xm(:),xpp(:),xpm(:),xmp(:),xmm(:),h(:)
    real(dp) :: f0
    integer :: i,j,n
    n=size(x)
    allocate(hessian(n,n),xp(n),xm(n),xpp(n),xpm(n),xmp(n),xmm(n),h(n))
    h=sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(x))
    f0=objective(x,context)
    do i=1,n
      xp=x;xm=x;xp(i)=xp(i)+h(i);xm(i)=xm(i)-h(i)
      hessian(i,i)=(objective(xp,context)-2.0_dp*f0+objective(xm,context))/(h(i)*h(i))
      do j=i+1,n
        xpp=x;xpm=x;xmp=x;xmm=x
        xpp(i)=xpp(i)+h(i);xpp(j)=xpp(j)+h(j)
        xpm(i)=xpm(i)+h(i);xpm(j)=xpm(j)-h(j)
        xmp(i)=xmp(i)-h(i);xmp(j)=xmp(j)+h(j)
        xmm(i)=xmm(i)-h(i);xmm(j)=xmm(j)-h(j)
        hessian(i,j)=(objective(xpp,context)-objective(xpm,context)- &
          objective(xmp,context)+objective(xmm,context))/(4.0_dp*h(i)*h(j))
        hessian(j,i)=hessian(i,j)
      end do
    end do
  end subroutine finite_hessian

  subroutine safe_standard_errors(hessian,se)
    real(dp),intent(in)::hessian(:,:)
    real(dp),allocatable,intent(out)::se(:)
    real(dp),allocatable::inv(:,:)
    type(tsc_error)::error
    integer::i,n
    n=size(hessian,1);allocate(se(n))
    call invert_matrix(hessian,inv,error)
    if(.not.error%ok())then
      se=ieee_value(0.0_dp,ieee_quiet_nan)
    else
      do i=1,n
        se(i)=sqrt(abs(inv(i,i)))
      end do
    end if
  end subroutine safe_standard_errors

  real(dp) function empirical_quantile_type7(x,p) result(q)
    real(dp),intent(in)::x(:),p
    real(dp),allocatable::y(:)
    real(dp)::h,g
    integer::j,n
    n=size(x)
    if(n==0)then
      q=ieee_value(0.0_dp,ieee_quiet_nan);return
    end if
    allocate(y(n));y=x;call sort_real(y)
    if(p<=0.0_dp)then;q=y(1);return
    else if(p>=1.0_dp)then;q=y(n);return
    end if
    h=1.0_dp+real(n-1,dp)*p
    j=int(floor(h));g=h-real(j,dp)
    if(j>=n)then;q=y(n)
    else;q=(1.0_dp-g)*y(j)+g*y(j+1)
    end if
  end function empirical_quantile_type7

  subroutine sort_real(x)
    real(dp),intent(inout)::x(:)
    real(dp)::key
    integer::i,j
    do i=2,size(x)
      key=x(i);j=i-1
      do while(j>=1)
        if(x(j)<=key)exit
        x(j+1)=x(j);j=j-1
      end do
      x(j+1)=key
    end do
  end subroutine sort_real

  real(dp) function log_sum_exp(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    m=maxval(x)
    if(.not.ieee_is_finite(m))then
      value=m
    else
      value=m+log(sum(exp(x-m)))
    end if
  end function log_sum_exp
end module tscopula_math
